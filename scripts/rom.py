"""
rom.py - Phase 5: thermisches Ersatznetzwerk (ROM) aus Fluent-Laeufen, Python-Teil.

Knoten: Umrichter, Netzteil, Schuetz, Trafo (Bauteil-Mitteltemperatur), Luft (Mittel), Wand (Mittel).
Waermestroeme (W):
  Bauteil -> Luft : q = gL_i * (dT^2 + 1 K^2)^(1/6) * dT (freie Konvektion, turbulent: h ~ dT^(1/3);
                    geglaettet, damit die Jacobi-Matrix bei dT = 0 (unbeheiztes Bauteil) regulaer bleibt)
  Bauteil -> Wand : q = gW_i * sigma * (T_i^4 - T_W^4)   (Strahlung, gW = eps_eff * A * F)
  Luft    -> Wand : q = gLW * (dT^2 + 1 K^2)^(1/6) * dT
  Wand -> Umgebung: wie die Fluent-Randbedingung: sum(h_s A_s) dT + eps sigma A (T_W^4 - T_U^4)
                    (h_s aus ergebnis.json des jeweiligen Laufs; Boden adiabat) -> NICHT angepasst
Identifikation: 9 Koeffizienten (gL, gW, gLW) per Least-Squares gegen 5 stationaere CFD-Laeufe
  (basis_grob + je ein Bauteil allein beheizt). Validierung: rom_val (andere Leistungen, 40 degC).
Tmax_i = T_i + r_i P_i  (r_i aus basis_grob: innerer Waermeleitwiderstand, linear).
Transient: C_i = rho c V (physikalisch, NICHT angepasst), C_Luft, C_Wand (2 mm Stahl).
           Vergleich mit dem transienten Fluent-Lauf (rom_transient).
Export: ergebnisse/rom/rom_parameter.json + linearisiertes Zustandsraummodell (A, B, C, D) fuer
        MATLAB (rom_matlab.m) und den Simulink-Generator (rom_simulink.m).

Aufruf: python rom.py
"""
import csv
import json
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import fsolve, least_squares

import config as C
import plotstil as PS

HIER = Path(__file__).resolve().parent
ZIEL = HIER / "ergebnisse" / "rom"
SIGMA = 5.670374e-8
BT = list(C.KOMPONENTEN)
KNOTEN = BT + ["Luft", "Wand"]
FIT = ["basis_grob"] + [f"rom_{n.lower()}" for n in BT]
VAL = ["rom_val"]
SEITE_A = {"decke": C.B * C.T, "links": C.H * C.T, "rechts": C.H * C.T, "rueck": C.B * C.H, "tuer": C.B * C.H}
A_AB = sum(SEITE_A.values())                                   # Flaeche ohne Boden
A_WAND = A_AB + C.B * C.T                                      # mit Boden (Waermekapazitaet)
V_LUFT = C.B * C.H * C.T - sum(np.prod(C.BLOCK_SIZE[n]) for n in BT)


def K(c):
    return c + 273.15


def lade(fall):
    return json.loads((HIER / "laeufe" / fall / "ergebnis.json").read_text(encoding="utf-8"))


class Netz:
    def __init__(self, theta):
        th = np.exp(theta[:9])
        self.gL, self.gW, self.gLW = th[0:4], th[4:8], th[8]
        # Schichtungsfaktoren: lokale Lufttemperatur am Bauteil = T_L + s_i (T_L - T_W)
        self.s = np.asarray(theta[9:13]) if len(theta) > 9 else np.zeros(4)

    def fluesse(self, T, P, T_U, h, eps_a=C.EMISSIVITAET):
        """Nettowaermebilanz je Knoten (W) bei Temperaturen T (K)."""
        Ti, TL, TW = T[:4], T[4], T[5]
        dL = Ti - (TL + self.s * (TL - TW))
        qiL = self.gL * (dL * dL + 1.0) ** (1 / 6) * dL
        qiW = self.gW * SIGMA * (Ti ** 4 - TW ** 4)
        dLW = TL - TW
        qLW = self.gLW * (dLW * dLW + 1.0) ** (1 / 6) * dLW
        hA = sum(h[s] * SEITE_A[s] for s in SEITE_A)
        qWU = hA * (TW - T_U) + eps_a * SIGMA * A_AB * (TW ** 4 - T_U ** 4)
        r = np.empty(6)
        r[:4] = P - qiL - qiW
        r[4] = qiL.sum() - qLW
        r[5] = qiW.sum() + qLW - qWU
        return r

    def stationaer(self, P, T_U, h):
        T0 = np.full(6, T_U + 15.0)
        T, info, ok, msg = fsolve(lambda T: self.fluesse(T, P, T_U, h), T0, full_output=True, xtol=1e-10)
        if ok != 1:
            raise RuntimeError(f"ROM stationaer nicht geloest: {msg}")
        return T


def daten(fall):
    e = lade(fall)
    p = e["parameter"]
    P = np.array([p[f"P_{n}"] for n in BT])
    T_cfd = np.array([e["Tmittel_C"][n] for n in BT] + [e["Tluft_mittel_C"], e["Tw_cfd_C"]])
    return dict(P=P, T_U=K(p["T_umg"]), h=e["h_aussen"], T=K(T_cfd), e=e)


def anpassen(faelle):
    d = [daten(f) for f in faelle]

    def res(theta):
        net = Netz(theta)
        out = []
        for x in d:
            try:
                out.append(net.stationaer(x["P"], x["T_U"], x["h"]) - x["T"])
            except RuntimeError:
                out.append(np.full(6, 50.0))
        return np.concatenate(out)

    # Hinweis: Ein Versuch mit 4 zusaetzlichen Schichtungsfaktoren (lokale Lufttemperatur je Bauteil)
    # verbesserte die Validierung nur von 1,14 auf 1,10 K bei physikalisch unplausiblen Werten
    # (Ueberanpassung) -> verworfen; s bleibt 0 (docs/log.md).
    theta0 = np.log(np.r_[np.full(4, 0.5), np.full(4, 0.1), 2.0])
    sol = least_squares(res, theta0, method="trf", x_scale=1.0, xtol=1e-12, ftol=1e-12)
    return Netz(sol.x), sol


def kapazitaeten():
    Cb = [C.BAUTEIL_MATERIAL[n]["rho"] * C.BAUTEIL_MATERIAL[n]["cp"] * np.prod(C.BLOCK_SIZE[n]) for n in BT]
    C_luft = 1.1 * 1006.0 * V_LUFT
    C_wand = C.WAND_RHO * C.WAND_CP * C.WAND * A_WAND
    return np.array(Cb + [C_luft, C_wand])


def linearisieren(net, P, T_U, h, Tss):
    """A, B fuer dx/dt = A x + B u mit x = T - Tss, u = [P_1..P_4, T_U] - u0 (numerische Ableitung)."""
    Cap = kapazitaeten()
    f0 = net.fluesse(Tss, P, T_U, h)
    A = np.zeros((6, 6))
    for j in range(6):
        dT = np.zeros(6)
        dT[j] = 1e-3
        A[:, j] = (net.fluesse(Tss + dT, P, T_U, h) - f0) / 1e-3 / Cap
    B = np.zeros((6, 5))
    for j in range(4):
        dP = np.zeros(4)
        dP[j] = 1e-3
        B[:, j] = (net.fluesse(Tss, P + dP, T_U, h) - f0) / 1e-3 / Cap
    B[:, 4] = (net.fluesse(Tss, P, T_U + 1e-3, h) - f0) / 1e-3 / Cap
    return A, B


def main():
    PS.stil()
    plt = PS.plt
    ZIEL.mkdir(parents=True, exist_ok=True)
    log = []
    t0 = time.time()
    net, sol = anpassen(FIT)
    t_fit = time.time() - t0
    basis = daten("basis_grob")
    e = basis["e"]
    r_int = {n: (e["Tmax_C"][n] - e["Tmittel_C"][n]) / e["parameter"][f"P_{n}"] for n in BT}
    log.append(f"ROM angepasst an {FIT} in {t_fit:.1f} s; RMS-Rest {np.sqrt(np.mean(sol.fun ** 2)):.3f} K")
    log.append("Koeffizienten: gL [W/K^4/3] " + ", ".join(f"{n} {g:.3f}" for n, g in zip(BT, net.gL))
               + " | gW [m2] " + ", ".join(f"{n} {g:.4f}" for n, g in zip(BT, net.gW))
               + f" | gLW {net.gLW:.3f} | Schichtung s " + ", ".join(f"{n} {v:+.3f}" for n, v in zip(BT, net.s)))

    # --- stationaerer Vergleich (Anpassung + Validierung)
    zeilen = []
    for fall in FIT + VAL:
        try:
            x = daten(fall)
        except FileNotFoundError:
            continue
        t1 = time.perf_counter()
        T = net.stationaer(x["P"], x["T_U"], x["h"])
        dt_rom = time.perf_counter() - t1
        for i, k in enumerate(KNOTEN):
            z = dict(fall=fall, rolle="Anpassung" if fall in FIT else "Validierung", knoten=k,
                     T_cfd_C=round(x["T"][i] - 273.15, 2), T_rom_C=round(T[i] - 273.15, 2),
                     fehler_K=round(T[i] - x["T"][i], 2))
            if i < 4:
                tmax_rom = T[i] - 273.15 + r_int[k] * x["P"][i]
                z.update(Tmax_cfd_C=round(x["e"]["Tmax_C"][k], 2), Tmax_rom_C=round(tmax_rom, 2),
                         Tmax_fehler_K=round(tmax_rom - x["e"]["Tmax_C"][k], 2))
            zeilen.append(z)
        fehler = [abs(T[i] - x["T"][i]) for i in range(6)]
        log.append(f"  {fall:14s} ({'Anpassung' if fall in FIT else 'VALIDIERUNG'}): max |Fehler| "
                   f"{max(fehler):.2f} K, ROM-Rechenzeit {1e3 * dt_rom:.1f} ms, CFD {x['e']['laufzeit_s'] / 60:.1f} min")
    felder = []
    for z in zeilen:
        felder += [k for k in z if k not in felder]
    with open(ZIEL / "rom_stationaer.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=felder)
        wr.writeheader()
        wr.writerows(zeilen)

    # --- Parameter + linearisiertes Modell exportieren
    Cap = kapazitaeten()
    Tss = net.stationaer(basis["P"], basis["T_U"], basis["h"])
    A, B = linearisieren(net, basis["P"], basis["T_U"], basis["h"], Tss)
    par = dict(knoten=KNOTEN, gL=net.gL.tolist(), gW=net.gW.tolist(), gLW=net.gLW, s=net.s.tolist(), r_int=r_int,
               C=Cap.tolist(), A_seite=SEITE_A, A_ab=A_AB, eps_aussen=C.EMISSIVITAET, sigma=SIGMA,
               h_basis=basis["h"], P_basis=basis["P"].tolist(), T_U_basis=basis["T_U"], T_ss_basis=Tss.tolist(),
               A_lin=A.tolist(), B_lin=B.tolist(), eingaenge=[f"P_{n}" for n in BT] + ["T_umg"])
    (ZIEL / "rom_parameter.json").write_text(json.dumps(par, indent=2), encoding="utf-8")
    log.append("Zeitkonstanten (linearisiert, Basis): "
               + ", ".join(f"{-1 / ev.real / 60:.1f} min" for ev in sorted(np.linalg.eigvals(A), key=lambda v: v.real)))

    # --- transient: ROM (nichtlinear) gegen Fluent
    tr = HIER / "laeufe" / "rom_transient" / "transient.csv"
    if tr.exists():
        import pandas as pd
        d = pd.read_csv(tr)
        h_tr = json.loads((HIER / "laeufe" / "rom_transient" / "ergebnis.json").read_text())["h_aussen"]
        T_U = basis["T_U"]
        P = basis["P"]
        t_end = d["t"].iloc[-1]
        t1 = time.perf_counter()
        ivp = solve_ivp(lambda t, T: net.fluesse(T, P, T_U, h_tr) / Cap, (0, t_end), np.full(6, T_U),
                        method="Radau", t_eval=d["t"].to_numpy(), rtol=1e-6, atol=1e-6)
        dt_rom = time.perf_counter() - t1
        cfd_lauf = json.loads((HIER / "laeufe" / "rom_transient" / "ergebnis.json").read_text())["laufzeit_s"]
        spalten = {n: f"tavg-{n.lower()}" for n in BT}
        spalten["Luft"] = "tavg-luft"
        fig, ax = plt.subplots(figsize=(7.0, 4.4))
        farben = dict(PS.BAUTEIL_FARBE, Luft=PS.NEUTRAL)
        fehl = {}
        for i, (k, c) in enumerate(spalten.items()):
            th = d["t"] / 3600
            ax.plot(th, d[c] - 273.15, color=farben[k], lw=2.0)
            ax.plot(ivp.t / 3600, ivp.y[i] - 273.15, color=farben[k], lw=1.4, ls="--")
            fehl[k] = float(np.max(np.abs(ivp.y[i] - d[c].to_numpy())))
            PS.endlabel(ax, th.iloc[-1], d[c].iloc[-1] - 273.15, k, farben[k])
        ax.plot([], [], color=PS.TEXT2, lw=2.0, label="Fluent (CFD)")
        ax.plot([], [], color=PS.TEXT2, lw=1.4, ls="--", label="ROM (RC-Netzwerk)")
        ax.legend(fontsize=8, loc="lower right")
        ax.set_xlim(0, t_end / 3600 * 1.15)
        ax.set_xlabel("Zeit in h")
        ax.set_ylabel("mittlere Temperatur in °C")
        ax.set_title(f"Aufheizen nach Einschalten (alle Verluste, 35 °C): CFD {cfd_lauf / 60:.0f} min "
                     f"Rechenzeit, ROM {1e3 * dt_rom:.0f} ms", fontsize=9)
        fig.savefig(ZIEL / "rom_transient.png")
        plt.close(fig)
        log.append("Transient max |ROM - CFD|: " + ", ".join(f"{k} {v:.2f} K" for k, v in fehl.items())
                   + f" | Rechenzeit ROM {1e3 * dt_rom:.0f} ms, CFD {cfd_lauf / 60:.1f} min")
        with open(ZIEL / "rom_transient.csv", "w", newline="", encoding="utf-8") as fh:
            wr = csv.writer(fh)
            wr.writerow(["t_s"] + [f"{k}_cfd_C" for k in spalten] + [f"{k}_rom_C" for k in spalten])
            for j, t in enumerate(ivp.t):
                wr.writerow([t] + [round(d[c].iloc[j] - 273.15, 3) for c in spalten.values()]
                            + [round(ivp.y[i][j] - 273.15, 3) for i in range(5)])
    # Bild stationaer
    import pandas as pd
    z = pd.DataFrame(zeilen)
    z = z[z.knoten.isin(BT)]
    fig, ax = plt.subplots(figsize=(5.2, 5.0))
    lo, hi = z[["T_cfd_C", "T_rom_C"]].min().min() - 2, z[["T_cfd_C", "T_rom_C"]].max().max() + 2
    ax.plot([lo, hi], [lo, hi], color=PS.GITTER, lw=1)
    for k in BT:
        q = z[z.knoten == k]
        a = q[q.rolle == "Anpassung"]
        v = q[q.rolle == "Validierung"]
        ax.plot(a.T_cfd_C, a.T_rom_C, "o", color=PS.BAUTEIL_FARBE[k], ms=7, label=k)
        ax.plot(v.T_cfd_C, v.T_rom_C, "D", color=PS.BAUTEIL_FARBE[k], ms=8, mfc=PS.FLAECHE, mew=2)
    ax.plot([], [], "D", color=PS.TEXT2, mfc=PS.FLAECHE, mew=2, label="Validierung (nicht angepasst)")
    ax.set_xlim(lo, hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("CFD: Bauteil-Mitteltemperatur in °C")
    ax.set_ylabel("ROM in °C")
    ax.set_title("ROM stationaer gegen CFD (Netz grob)", fontsize=9)
    ax.legend(fontsize=8)
    fig.savefig(ZIEL / "rom_stationaer.png")
    plt.close(fig)
    (HIER / "logs" / "rom.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))


if __name__ == "__main__":
    main()
