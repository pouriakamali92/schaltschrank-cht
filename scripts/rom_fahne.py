"""
rom_fahne.py - ROM-Erweiterung: Fahnen-/Schichtungskopplung zwischen den Bauteilen (Einflusskoeffizienten).

Basis-ROM (rom.py): jedes Bauteil tauscht Waerme mit der MITTLEREN Lufttemperatur aus. Erweiterung:
    T_lokal_i = T_L + sum_j K_ij P_j      (K_ii = 0, K_ij in K/W; lineare Ueberlagerung = Annahme)
    qL_i      = gL_i ((T_i - T_lokal_i)^2 + 1)^(1/6) (T_i - T_lokal_i)
Die Waerme geht weiter in den Luftknoten -> Energiebilanz bleibt exakt.
Varianten: A ohne Fahne (= rom.py), B volle 4x4-Matrix (12 K, Vorzeichen frei, mit Grundkoeff.),
           C nur geometrisch plausible Fahnenpaare (Bauteil oberhalb, x-Ueberlappung), K >= 0, Grundkoeff. fest,
           D wie C, aber gemeinsam mit den 9 Grundkoeffizienten angepasst.
Bewertung (gegen Ueberanpassung): Leave-one-out ueber die 5 Anpassungslaeufe, unabhaengige Validierung rom_val,
           transienter Lauf. Uebernahme nur, wenn Validierung UND Leave-one-out um > 0,1 K besser als A.
Ausgabe: ergebnisse/rom/rom_fahne.csv, rom_fahne_K.csv, rom_fahne.png, logs/rom_fahne.txt

Aufruf: python rom_fahne.py
"""
import csv
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import solve_ivp
from scipy.optimize import least_squares

import config as C
import plotstil as PS
import rom as R

HIER = Path(__file__).resolve().parent
ZIEL = HIER / "ergebnisse" / "rom"
BT = R.BT
AUS = ~np.eye(4, dtype=bool)               # Positionen der 12 Nebendiagonal-Koeffizienten


class NetzFahne(R.Netz):
    def __init__(self, theta):
        super().__init__(theta[:9])
        self.K = np.zeros((4, 4))
        if len(theta) > 9:
            self.K[AUS] = theta[9:21]

    def fluesse(self, T, P, T_U, h, eps_a=C.EMISSIVITAET):
        Ti, TL, TW = T[:4], T[4], T[5]
        dL = Ti - (TL + self.K @ P)
        qiL = self.gL * (dL * dL + 1.0) ** (1 / 6) * dL
        qiW = self.gW * R.SIGMA * (Ti ** 4 - TW ** 4)
        dLW = TL - TW
        qLW = self.gLW * (dLW * dLW + 1.0) ** (1 / 6) * dLW
        hA = sum(h[s] * R.SEITE_A[s] for s in R.SEITE_A)
        qWU = hA * (TW - T_U) + eps_a * R.SIGMA * R.A_AB * (TW ** 4 - T_U ** 4)
        r = np.empty(6)
        r[:4] = P - qiL - qiW
        r[4] = qiL.sum() - qLW
        r[5] = qiW.sum() + qLW - qWU
        return r


def fahnenpaare(rand=0.05):
    """A priori aus der Geometrie: Fahne von j wirkt auf i, wenn i oberhalb von j liegt (Unterkante i >= Oberkante j)
    und sich die x-Bereiche bis auf `rand` ueberlappen. Liefert Indizes in theta[9:21]."""
    paare, idx = [], []
    pos = [(i, j) for i in range(4) for j in range(4) if i != j]          # Reihenfolge wie K[AUS]
    for k, (i, j) in enumerate(pos):
        oi, si, oj, sj = C.BLOCK_ORIGIN[BT[i]], C.BLOCK_SIZE[BT[i]], C.BLOCK_ORIGIN[BT[j]], C.BLOCK_SIZE[BT[j]]
        oben = oi[1] >= oj[1] + sj[1] - 1e-9
        x_ueberlapp = oi[0] < oj[0] + sj[0] + rand and oj[0] < oi[0] + si[0] + rand
        if oben and x_ueberlapp:
            paare.append(f"{BT[j]}->{BT[i]}")
            idx.append(9 + k)
    return paare, idx


def anpassen(faelle, frei, theta_ref, untere=None):
    """Passt die Eintraege `frei` von theta (9 log-Grundkoeff. + 12 K) an, Rest fest auf theta_ref."""
    d = [R.daten(f) for f in faelle]
    frei = np.asarray(frei)

    def voll(p):
        th = theta_ref.copy()
        th[frei] = p
        return th

    def res(p):
        net = NetzFahne(voll(p))
        out = []
        for x in d:
            try:
                out.append(net.stationaer(x["P"], x["T_U"], x["h"]) - x["T"])
            except RuntimeError:
                out.append(np.full(6, 50.0))
        return np.concatenate(out)

    lo = np.full(len(frei), -np.inf) if untere is None else untere
    p0 = np.clip(theta_ref[frei], lo, np.inf)
    sol = least_squares(res, p0, bounds=(lo, np.inf), method="trf", xtol=1e-12, ftol=1e-12)
    return voll(sol.x), sol


def fehler(net, fall):
    x = R.daten(fall)
    return net.stationaer(x["P"], x["T_U"], x["h"]) - x["T"]


def transient(net):
    d = pd.read_csv(HIER / "laeufe" / "rom_transient" / "transient.csv")
    h = json.loads((HIER / "laeufe" / "rom_transient" / "ergebnis.json").read_text())["h_aussen"]
    b = R.daten("basis_grob")
    Cap = R.kapazitaeten()
    iv = solve_ivp(lambda t, T: net.fluesse(T, b["P"], b["T_U"], h) / Cap, (0, d["t"].iloc[-1]),
                   np.full(6, b["T_U"]), method="Radau", t_eval=d["t"].to_numpy(), rtol=1e-8, atol=1e-8)
    sp = ["tavg-umrichter", "tavg-netzteil", "tavg-schuetz", "tavg-trafo", "tavg-luft"]
    return {k: float(np.max(np.abs(iv.y[i] - d[c].to_numpy()))) for i, (k, c) in enumerate(zip(BT + ["Luft"], sp))}, iv, d


def bewerte(name, frei, untere, theta_ref):
    """Anpassung, Validierung, Leave-one-out und transienter Vergleich einer Variante."""
    th, sol = anpassen(R.FIT, frei, theta_ref, untere)
    net = NetzFahne(th)
    e = {f: fehler(net, f) for f in R.FIT + R.VAL}
    loo = {}
    for weg in R.FIT:
        tl, _ = anpassen([f for f in R.FIT if f != weg], frei, th, untere)
        loo[weg] = float(np.max(np.abs(fehler(NetzFahne(tl), weg))))
    tr, iv, _ = transient(net)
    return dict(name=name, n=len(frei), th=th, net=net, rms=float(np.sqrt(np.mean(sol.fun ** 2))), e=e, loo=loo,
                tr=tr, iv=iv,
                fit=max(float(np.max(np.abs(e[f]))) for f in R.FIT),
                val=float(np.max(np.abs(e["rom_val"]))), loo_max=max(loo.values()), loo_mittel=float(np.mean(list(loo.values()))),
                tr_max=max(tr.values()))


def main():
    PS.stil()
    plt = PS.plt
    log = []
    t0 = time.time()
    # Referenz = Basis-ROM aus rom.py (9 Grundkoeffizienten, K = 0)
    th0, _ = anpassen(R.FIT, list(range(9)), np.r_[np.log(np.r_[np.full(4, 0.5), np.full(4, 0.1), 2.0]), np.zeros(12)])
    paare, idx3 = fahnenpaare()
    nullK = np.zeros(len(idx3))
    var = [
        bewerte("A ohne Fahne (rom.py)", list(range(9)), None, th0),
        bewerte("B volle Matrix (9+12, Vorzeichen frei)", list(range(21)), None, th0),
        bewerte(f"C {len(idx3)} Fahnenpaare, K>=0, Grundkoeff. fest", idx3, nullK, th0),
        bewerte(f"D {len(idx3)} Fahnenpaare, K>=0, gemeinsam (9+{len(idx3)})", list(range(9)) + idx3,
                np.r_[np.full(9, -np.inf), nullK], th0),
    ]
    log.append(f"ROM-Fahnenkopplung ({time.strftime('%Y-%m-%d %H:%M')}, {time.time() - t0:.0f} s). "
               f"Geometrische Fahnenpaare (Bauteil oberhalb, x-Ueberlappung +-5 cm): {', '.join(paare)}")
    log.append("Variante                                       Koeff  RMS-Anp.  max Anp.  LOO max/mittel  VALID. rom_val  transient")
    for v in var:
        log.append(f"  {v['name']:44s} {v['n']:4d}  {v['rms']:6.3f} K  {v['fit']:6.2f} K  {v['loo_max']:5.2f}/{v['loo_mittel']:4.2f} K"
                   f"     {v['val']:5.2f} K      {v['tr_max']:5.2f} K")
    for v in var[1:]:
        K = v["net"].K
        log.append(f"  K [K/W] {v['name'][:1]}: " + ", ".join(f"{BT[j]}->{BT[i]} {K[i, j]:+.4f}"
                                                         for i in range(4) for j in range(4) if i != j and abs(K[i, j]) > 1e-6))
    for v in var:
        log.append(f"  rom_val {v['name'][:1]}: " + " ".join(f"{k} {x:+.2f}" for k, x in zip(R.KNOTEN, v["e"]["rom_val"])))

    # Entscheidung nach festem Kriterium: nur uebernehmen, wenn Validierung UND Leave-one-out besser als A
    a = var[0]
    besser = [v for v in var[1:] if v["val"] < a["val"] - 0.1 and v["loo_max"] < a["loo_max"] - 0.1
              and v["tr_max"] <= a["tr_max"] + 0.1]
    if besser:
        w = min(besser, key=lambda v: v["loo_max"])
        log.append(f"ERGEBNIS: Variante {w['name']} verbessert Validierung ({a['val']:.2f} -> {w['val']:.2f} K) und "
                   f"Leave-one-out ({a['loo_max']:.2f} -> {w['loo_max']:.2f} K).")
    else:
        log.append("ERGEBNIS: Keine Variante verbessert Validierung UND Leave-one-out gegenueber dem Basis-ROM -> "
                   "Fahnenkopplung NICHT uebernommen; Basis-ROM (rom.py, rom_export) bleibt gueltig. Mit 5 stationaeren "
                   "Laeufen (davon 4 Einzelbeheizungen) ist die Kopplung nicht belastbar identifizierbar.")

    with open(ZIEL / "rom_fahne.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["variante", "koeffizienten", "rms_anpassung_K", "max_anpassung_K", "loo_max_K", "loo_mittel_K",
                     "validierung_rom_val_K", "transient_max_K"] + [f"loo_{f}_K" for f in R.FIT])
        for v in var:
            wr.writerow([v["name"], v["n"], round(v["rms"], 3), round(v["fit"], 3), round(v["loo_max"], 3),
                         round(v["loo_mittel"], 3), round(v["val"], 3), round(v["tr_max"], 3)]
                        + [round(v["loo"][f], 3) for f in R.FIT])
    with open(ZIEL / "rom_fahne_K.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["variante", "beheizt_j", "betroffen_i", "K_ij_K_per_W"])
        for v in var[1:]:
            for i in range(4):
                for j in range(4):
                    if i != j and abs(v["net"].K[i, j]) > 1e-9:
                        wr.writerow([v["name"], BT[j], BT[i], round(v["net"].K[i, j], 5)])

    # Bild: drei Fehlermasse je Variante (gruppierte Balken, eine Achse)
    fig, ax = plt.subplots(figsize=(7.4, 3.8))
    masse = [("Anpassung (max)", "fit", PS.NEUTRAL), ("Leave-one-out (max)", "loo_max", PS.KATEGORIE[1]),
             ("Validierung rom_val", "val", PS.KATEGORIE[0]), ("transient 8 h", "tr_max", PS.KATEGORIE[2])]
    x = np.arange(len(var))
    for k, (lab, key, col) in enumerate(masse):
        ax.bar(x + (k - 1.5) * 0.2, [v[key] for v in var], 0.18, color=col, label=lab)
    ax.set_xticks(x)
    ax.set_xticklabels([v["name"].replace(" (", "\n(").replace(", ", "\n", 1) for v in var], fontsize=7)
    ax.set_ylabel("max. |ROM − CFD| in K")
    ax.set_title("ROM-Fahnenkopplung: " + ("Vorhersage verbessert" if besser else
                 "bessere Anpassung, aber keine bessere Vorhersage"), fontsize=9)
    ax.legend(fontsize=7, ncol=2)
    fig.savefig(ZIEL / "rom_fahne.png")
    plt.close(fig)
    (HIER / "logs" / "rom_fahne.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))


if __name__ == "__main__":
    main()


