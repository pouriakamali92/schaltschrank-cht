"""
plausibilitaet_audit.py - unabhaengige Nachpruefung der Luecken (aendert KEINE Ergebnisse).

Liest vorhandene Ergebnisse (laeufe/*/ergebnis.json, logs/gci.txt, ergebnisse/*, rom_export/) und die nur
gelesenen Fluent-Groessen aus logs/audit_fluent.json (erzeugt von audit_fluent_lesen.py: Bauteil-
Oberflaechentemperatur, Strahlungs-/Gesamtwaermestrom je Bauteil, Luefter-Massenstrom und Austrittstemperatur).
Schreibt nur logs/plausibilitaet_audit.txt.

  1 Innentemperaturen: Handmodell (Churchill-Chu / horizontale Platten + Strahlung zur Wand, Wand aus
    handrechnung.py) mit Bandbreite -> Oberflaechen-, Luft-, Tmax-Temperaturen gegen CFD
  2 GCI nach Celik (2008) aus den Rohwerten selbst nachgerechnet
  3 EM: Skintiefe, R_ac/R_dc >= 1 und monoton, Kopplung dP ~ alpha dT
  4 ROM: C = m cp, Anfangssteigung P/C gegen CFD, tau = C/G, 63-%-Zeiten, Steady State gegen CFD
  5 Luefter: P_Luft = mdot cp dT gegen Bilanzrest

Aufruf: python plausibilitaet_audit.py
"""
import itertools
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import fsolve

import config as C
import em_modell as EM

HIER = Path(__file__).resolve().parent
SIG = 5.670374e-8
OUT = []
ZAEHLER = {"OK": 0, "FRAGLICH": 0, "FEHLER": 0}


def zeile(was, wert, quelle, erwartung, urteil, notiz=""):
    ZAEHLER[urteil] += 1
    OUT.append(f"  [{urteil:8s}] {was}: {wert} | Quelle: {quelle} | Erwartung: {erwartung}" + (f" | {notiz}" if notiz else ""))


def text(s=""):
    OUT.append(s)


def lade(fall):
    return json.loads((HIER / "laeufe" / fall / "ergebnis.json").read_text(encoding="utf-8"))


# =============================================================================================
# 1 Innentemperaturen
# =============================================================================================
def nu_vertikal(Ra, Pr=C.LUFT_PR):
    return (0.825 + 0.387 * Ra ** (1 / 6) / (1 + (0.492 / Pr) ** (9 / 16)) ** (8 / 27)) ** 2


def nu_oben_warm(Ra):            # Oberseite warme Platte / Unterseite kalte Platte
    return 0.54 * Ra ** 0.25 if Ra < 1e7 else 0.15 * Ra ** (1 / 3)


def nu_unten_warm(Ra):           # Unterseite warme Platte / Oberseite kalte Platte
    return 0.52 * Ra ** 0.2


def ra(dT, L, Tf):
    return C.G * (1 / Tf) * abs(dT) * L ** 3 / (C.LUFT_NU * C.LUFT_ALPHA)


def h_konv(dT, L, art, Tf, fc):
    if abs(dT) < 1e-6:
        return 0.0
    R = ra(dT, L, Tf)
    nu = {"v": nu_vertikal, "o": nu_oben_warm, "u": nu_unten_warm}[art](R)
    return fc * nu * C.LUFT_LAMBDA / L


def flaechen(n):
    x, y, z = C.BLOCK_SIZE[n]
    return dict(front=(x * y, "v", y), back=(x * y, "v", y), seiten=(2 * z * y, "v", y),
                oben=(x * z, "o", x * z / (2 * (x + z))), unten=(x * z, "u", x * z / (2 * (x + z))))


def handmodell(Tw, fc, back, eps, F, fin):
    """Unbekannte: 4 Oberflaechentemperaturen + Luft. Wand gleichmaessig Tw (aus handrechnung.py)."""
    BT = list(C.KOMPONENTEN)
    A_in = 2 * C.B * C.H + 2 * C.T * C.H     # senkrechte Innenwaende
    A_de = A_bo = C.B * C.T

    def res(x):
        Ts, Ta = x[:4], x[4]
        r = np.empty(5)
        qk_sum = 0.0
        for i, n in enumerate(BT):
            qk = 0.0
            for name, (A, art, L) in flaechen(n).items():
                f = back if name == "back" else 1.0
                qk += f * h_konv(Ts[i] - Ta, L, art, (Ts[i] + Ta) / 2, fc) * A * (Ts[i] - Ta)
            A_ges = sum(a for a, _, _ in flaechen(n).values())
            qr = eps * F * SIG * A_ges * (Ts[i] ** 4 - Tw ** 4)
            r[i] = C.KOMPONENTEN[n] - qk - qr
            qk_sum += qk
        Tf = (Ta + Tw) / 2
        q_in = (h_konv(Ta - Tw, C.H, "v", Tf, fin) * A_in + h_konv(Ta - Tw, A_de / (2 * (C.B + C.T)), "o", Tf, fin) * A_de
                + h_konv(Ta - Tw, A_bo / (2 * (C.B + C.T)), "u", Tf, fin) * A_bo) * (Ta - Tw)
        r[4] = qk_sum - q_in
        return r

    x = fsolve(res, np.r_[np.full(4, Tw + 40), Tw + 6], xtol=1e-10)
    Ts, Ta = x[:4], x[4]
    anteil_rad = []
    for i, n in enumerate(BT):
        A_ges = sum(a for a, _, _ in flaechen(n).values())
        anteil_rad.append(eps * F * SIG * A_ges * (Ts[i] ** 4 - Tw ** 4) / C.KOMPONENTEN[n])
    return Ts, Ta, np.array(anteil_rad)


def box_dT(n):
    """Exakte Reihenloesung Quader mit gleichmaessiger Quelle, Oberflaeche isotherm: (Zentrum, Volumenmittel) - Oberflaeche."""
    a, b, c = C.BLOCK_SIZE[n]
    q = C.KOMPONENTEN[n] / (a * b * c)
    k = C.BAUTEIL_MATERIAL[n]["k"]
    zen = mit = 0.0
    for l in range(1, 60, 2):
        for m in range(1, 60, 2):
            for nn in range(1, 60, 2):
                lam = math.pi ** 2 * (l * l / a ** 2 + m * m / b ** 2 + nn * nn / c ** 2)
                koef = 64 / (math.pi ** 3 * l * m * nn * lam)
                zen += koef * (-1) ** ((l + m + nn - 3) // 2)
                mit += koef * 8 / (math.pi ** 3 * l * m * nn)
    return q / k * zen, q / k * mit


def teil1(fl):
    text("1) INNENTEMPERATUREN - unabhaengiges Handmodell (Basis 35 degC, Netz mittel)")
    text("   Modell: Bauteiloberflaeche -> Luft (Churchill-Chu senkrecht, 0,54/0,15 Ra^n oben, 0,52 Ra^0,2 unten),")
    text("   Bauteil -> Wand Strahlung eps_eff*F*sigma*(Ts^4-Tw^4); Luft -> Innenwand (gleiche Korrelationen);")
    text("   Wand Tw = 43,49 degC aus handrechnung.py (NICHT aus CFD). Innerer Anstieg: exakte Quader-Reihenloesung.")
    text("   Bandbreite (32 Ecken): Korrelation x0,8/x1,2; Rueckseite (2-cm-Spalt) Konvektion 0/voll; eps_eff 0,8/0,9;")
    text("   Sichtanteil zur Wand F 0,8/1,0; innen Luft->Wand x0,8/x1,2.")
    e = lade("basis_mittel")
    Tw = e["Tw_hand_C"] + 273.15
    BT = list(C.KOMPONENTEN)
    erg = []
    for fc, back, eps, F, fin in itertools.product([0.8, 1.2], [0.0, 1.0], [0.8, 0.9], [0.8, 1.0], [0.8, 1.2]):
        erg.append(handmodell(Tw, fc, back, eps, F, fin))
    nom = handmodell(Tw, 1.0, 0.5, 0.85, 0.9, 1.0)
    Ts_all = np.array([x[0] for x in erg]) - 273.15
    Ta_all = np.array([x[1] for x in erg]) - 273.15
    rad_all = np.array([x[2] for x in erg])
    ta_lo, ta_hi = Ta_all.min(), Ta_all.max()
    q = "Handmodell (plausibilitaet_audit.py) / CFD laeufe/basis_mittel/ergebnis.json"
    ok = ta_lo - 1 <= e["Tluft_mittel_C"] <= ta_hi + 1
    zeile("Luft-Mittel", f"Hand {nom[1] - 273.15:.1f} degC [{ta_lo:.1f}..{ta_hi:.1f}], CFD {e['Tluft_mittel_C']:.1f}", q,
          "CFD im Band (+-1 K)", "OK" if ok else "FRAGLICH",
          "Handmodell kennt keine Schichtung (CFD oben-unten 13 K)")
    b = fl.get("basis_mittel", {}) if fl else {}
    for i, n in enumerate(BT):
        lo, hi = Ts_all[:, i].min(), Ts_all[:, i].max()
        nk = n.lower()
        zen, mit = box_dT(n)
        if nk in b:
            Ts_cfd = b[nk][f"t-{nk}"] - 273.15
            urteil = "OK" if lo - 2 <= Ts_cfd <= hi + 2 else "FRAGLICH"
            zeile(f"{n} Oberflaeche", f"Hand {nom[0][i] - 273.15:.1f} degC [{lo:.1f}..{hi:.1f}], CFD {Ts_cfd:.1f} (flaechengemittelt)",
                  "Handmodell / logs/audit_fluent.json (aus basis_mittel.dat.h5 gelesen)", "CFD im Band (+-2 K)", urteil)
            qr, qt = b[nk][f"qr-{nk}"], b[nk][f"qt-{nk}(without-sources)"]
            ant = abs(qr) / abs(qt) if qt else float("nan")
            r_lo, r_hi = rad_all[:, i].min(), rad_all[:, i].max()
            zeile(f"{n} Strahlungsanteil", f"Hand {nom[2][i]:.0%} [{r_lo:.0%}..{r_hi:.0%}], CFD {ant:.0%} (|q_rad|/|q_ges| = {abs(qr):.1f}/{abs(qt):.1f} W)",
                  "Handmodell / logs/audit_fluent.json", "CFD im Band (+-10 %-Punkte)",
                  "OK" if r_lo - 0.10 <= ant <= r_hi + 0.10 else "FRAGLICH",
                  "CFD zaehlt auch Strahlungsaustausch mit anderen Bauteilen")
            for was, d_h, cfd in (("Tmax - Oberflaeche", zen, e["Tmax_C"][n] - Ts_cfd),
                                  ("Tmittel - Oberflaeche", mit, e["Tmittel_C"][n] - Ts_cfd)):
                urteil = "OK" if 0.7 * d_h <= cfd <= 1.3 * d_h else "FRAGLICH"
                zeile(f"{n} {was}", f"Reihenloesung {d_h:.2f} K, CFD {cfd:.2f} K ({cfd / d_h - 1:+.0%}, {cfd - d_h:+.2f} K)",
                      "Reihenloesung / ergebnis.json + audit_fluent.json", "+-30 %", urteil,
                      "CFD systematisch hoeher: Oberflaeche nicht isotherm (Rueckseite im 2-cm-Spalt waermer als Mittel)")
        tmax_lo, tmax_hi = lo + 0.7 * zen, hi + 1.3 * zen
        t = e["Tmax_C"][n]
        zeile(f"{n} Tmax", f"Hand [{tmax_lo:.1f}..{tmax_hi:.1f}] degC (nominal {nom[0][i] - 273.15 + zen:.1f}), CFD {t:.1f}",
              q, "CFD im Band", "OK" if tmax_lo <= t <= tmax_hi else "FRAGLICH",
              f"Bandbreite {tmax_hi - tmax_lo:.0f} K")
    text("   NICHT validiert damit: Stroemungsbild/Fahnen/Schichtung und lokale Hotspots; ob die Ersatzbloecke")
    text("   (homogenes k, gleichmaessige Quelle) reale Bauteile abbilden; die Verlustleistungen selbst; der")
    text("   Aussenwaermeuebergang (Handrechnung und CFD nutzen DIESELBEN Korrelationen). Es ist ein")
    text("   Groessenordnungs-/Konsistenztest zweier Modelle, keine Validierung gegen Messung. Die Tmax-Baender sind")
    text("   12-23 K breit: der Test erkennt grobe Fehler (Faktor, Vorzeichen, falsche Quelle), nicht 2-3 K.")
    text()


# =============================================================================================
# 2 GCI
# =============================================================================================
def celik(f1, f2, f3, r21, r32):
    e21, e32 = f2 - f1, f3 - f2
    if e21 == 0 or e32 == 0:
        return dict(p=float("nan"), mono=False, gci=float("nan"), s=0)
    s = np.sign(e32 / e21)
    p = abs(math.log(abs(e32 / e21))) / math.log(r21)
    for _ in range(200):
        try:
            qp = math.log((r21 ** p - s) / (r32 ** p - s))
        except ValueError:
            qp = 0.0
        p_neu = abs(math.log(abs(e32 / e21)) + qp) / math.log(r21)
        if abs(p_neu - p) < 1e-10:
            break
        p = p_neu
    ea21 = abs((f1 - f2) / f1)
    gci = 1.25 * ea21 / (r21 ** p - 1) if p > 0 else float("nan")
    return dict(p=p, mono=s > 0, gci=gci, s=s)


def teil2():
    text("2) GCI (Celik et al. 2008) aus den Rohwerten nachgerechnet")
    V = C.B * C.H * C.T
    N = {"grob": 45288, "mittel": 134913, "fein": 353186}
    h = {k: (V / v) ** (1 / 3) for k, v in N.items()}
    r21, r32 = h["mittel"] / h["fein"], h["grob"] / h["mittel"]
    zeile("Netzweite h fein/mittel/grob", f"{1e3 * h['fein']:.2f}/{1e3 * h['mittel']:.2f}/{1e3 * h['grob']:.2f} mm",
          "logs/gate1_netz.txt (Zellzahlen), V = 0,384 m3", "wie logs/gci.txt 10,3/14,2/20,4", "OK")
    zeile("r21 / r32", f"{r21:.3f} / {r32:.3f}", "berechnet", ">= 1,3 (Celik)",
          "OK" if min(r21, r32) >= 1.3 else "FRAGLICH")
    roh = {  # phi = Anstieg ueber Umgebung (K): grob, mittel, fein  [logs/gci.txt]
        "Tmax Umrichter": (43.034, 43.046, 43.043), "Tmax Netzteil": (52.646, 52.149, 52.359),
        "Tmax Schuetz": (69.247, 68.813, 69.558), "Tmax Trafo": (41.419, 41.308, 41.245),
        "Luft Mittel": (14.965, 15.039, 14.974), "Wand (Q-gew.)": (8.503, 8.501, 8.513)}
    for name, (g, m, f) in roh.items():
        c = celik(f, m, g, r21, r32)
        plausibel = 0.5 <= c["p"] <= 3
        urteil = "OK" if (c["mono"] and plausibel) else "FRAGLICH"
        grund = []
        if not c["mono"]:
            grund.append("nicht monoton (oszillierend)")
        if not plausibel:
            grund.append(f"p ausserhalb 0,5..3")
        zeile(name, f"phi g/m/f {g}/{m}/{f} K -> p {c['p']:.2f}, {'monoton' if c['mono'] else 'oszillierend'}, "
                    f"GCI_fein {100 * c['gci']:.2f} % ({c['gci'] * f:.3f} K), Spannweite {max(g, m, f) - min(g, m, f):.3f} K",
              "Rohwerte logs/gci.txt", "monoton und 0,5 <= p <= 3", urteil, "; ".join(grund))
    text("   Folge: Nur Trafo erfuellt die GCI-Voraussetzungen. Fuer alle anderen Groessen ist der GCI-Wert nicht")
    text("   belastbar; im Projekt wird deshalb die Spannweite (max. 0,75 K) als Netzunsicherheit angegeben.")
    text("   Die Spannweiten liegen in der Groesse der zeitlichen Schwankung (+-0,17..0,47 K, logs/gci.txt).")
    text()


# =============================================================================================
# 3 EM
# =============================================================================================
def teil3():
    text("3) EM - Skintiefe, R_ac/R_dc, Kopplung")
    rho20, alpha, mu0 = 1.7241e-8, 0.00393, 4e-7 * math.pi      # Cu weich: 1/58 Ohm mm2/m, alpha 0,00393 1/K
    for T in (20.0, 80.0):
        rho = rho20 * (1 + alpha * (T - 20))
        for f in (50.0, 1000.0, 10000.0):
            d_hand = math.sqrt(rho / (math.pi * f * mu0))
            d_code = EM.skintiefe(f, T)
            abw = d_code / d_hand - 1
            zeile(f"Skintiefe {f:g} Hz, {T:g} degC", f"Hand {1e3 * d_hand:.3f} mm, Code {1e3 * d_code:.3f} mm ({100 * abw:+.2f} %)",
                  "em_modell.skintiefe / Hand (rho20 1,7241e-8, alpha 0,00393)", "|Abw.| < 0,5 %",
                  "OK" if abs(abw) < 0.005 else "FEHLER")
    d50 = math.sqrt(rho20 * (1 + alpha * 60) / (math.pi * 50 * mu0))
    zeile("Bezug Gate 4: Skintiefe 50 Hz / 80 degC", f"{1e3 * d50:.1f} mm", "logs/gate4_em.txt nennt 10,4 mm", "10,4 mm",
          "OK" if abs(1e3 * d50 - 10.4) < 0.1 else "FEHLER")
    s = pd.read_csv(HIER / "ergebnisse" / "em" / "schienen_rac_rdc.csv")
    for c in s.columns[1:]:
        mn, mono = s[c].min(), bool(np.all(np.diff(s[c].to_numpy()) >= -1e-9))
        zeile(f"Schienen R_ac/R_dc {c}", f"min {mn:.6f}, monoton steigend: {mono}, max {s[c].max():.3f} bei {s.f_Hz.iloc[-1]:.0f} Hz",
              "ergebnisse/em/schienen_rac_rdc.csv", ">= 1 und monoton in f", "OK" if mn >= 1 - 1e-9 and mono else "FEHLER")
    d = pd.read_csv(HIER / "ergebnisse" / "em" / "trafo_dowell.csv")
    for c in d.columns[1:]:
        mn, mono = d[c].min(), bool(np.all(np.diff(d[c].to_numpy()) >= -1e-9))
        zeile(f"Dowell F_R {c} Lagen", f"min {mn:.6f}, monoton: {mono}", "ergebnisse/em/trafo_dowell.csv",
              ">= 1 und monoton in f", "OK" if mn >= 1 - 1e-9 and mono else "FEHLER")
    # Kopplung: Werte aus logs/em_kopplung.txt (Kopplung 1 bei 80 degC, Kopplung 3 = Endwerte)
    start = dict(L1=(80.00, 6.18), L2=(80.00, 6.19), L3=(80.00, 6.18))
    ende = dict(L1=(65.09, 5.90), L2=(66.64, 5.94), L3=(65.15, 5.90))
    for k in start:
        (T0, P0), (T1, P1) = start[k], ende[k]
        erw = (1 + alpha * (T1 - 20)) / (1 + alpha * (T0 - 20))
        ist = P1 / P0
        zeile(f"Kopplung Schiene {k}", f"T {T0:.1f} -> {T1:.1f} degC, P {P0:.2f} -> {P1:.2f} W, P1/P0 {ist:.4f}; "
                                    f"nur R_dc(T): {erw:.4f}", "logs/em_kopplung.txt",
              "P sinkt, weil T unter dem Startwert 80 degC endet; P1/P0 ~ R(T1)/R(T0) (Abw. < 1 %, Rest: AC-Faktor/Rundung)",
              "OK" if P1 < P0 and abs(ist / erw - 1) < 0.01 else "FRAGLICH")
    text("   Hinweis zur Vorgabe 'Endverluste >= Anfangsverluste': gilt nur, wenn die Startschaetzung UNTER der")
    text("   Endtemperatur liegt. Hier Start 80 degC, Ende ~65 degC -> Verluste muessen sinken (beobachtet -4,5 %).")
    zeile("Kopplung Wicklung", "T 80,0 -> 79,4 degC, P 59,49 -> 59,50 W", "logs/em_kopplung.txt",
          "nahezu konstant (starker Skineffekt: P ~ sqrt(rho), dazu dT nur 0,6 K)", "OK")
    text()


# =============================================================================================
# 4 ROM
# =============================================================================================
def teil4():
    text("4) ROM - Kapazitaeten, Zeitkonstanten, Steady State")
    par = json.loads((HIER / "rom_export" / "rom_parameter.json").read_text(encoding="utf-8"))
    BT = list(C.KOMPONENTEN)
    Cap = par["C"]["wert"]
    G = np.array(par["G_matrix"]["wert"])
    tr = pd.read_csv(HIER / "laeufe" / "rom_transient" / "transient.csv")
    rt = pd.read_csv(HIER / "ergebnisse" / "rom" / "rom_transient.csv")
    for i, n in enumerate(BT):
        m = C.BAUTEIL_MATERIAL[n]["rho"] * np.prod(C.BLOCK_SIZE[n])
        c_hand = m * C.BAUTEIL_MATERIAL[n]["cp"]
        zeile(f"{n} C = m cp", f"m {m:.2f} kg, C {c_hand:.0f} J/K, ROM {Cap[i]:.0f} J/K", "config.py / rom_export/rom_parameter.json",
              "gleich", "OK" if abs(c_hand / Cap[i] - 1) < 1e-9 else "FEHLER")
        # Anfangssteigung: in den ersten 300 s fast adiabat -> dT ~ P t / C (minus Abgabe)
        col = f"tavg-{n.lower()}"
        dT_cfd = tr[col].iloc[1] - tr[col].iloc[0]
        t1 = tr["t"].iloc[1]
        Gi = G[i, 4] + G[i, 5]
        tau = Cap[i] / Gi
        dT_erw = C.KOMPONENTEN[n] / Gi * (1 - math.exp(-t1 / tau))
        zeile(f"{n} Anstieg erste {t1:.0f} s", f"CFD {dT_cfd:.2f} K, adiabat P t/C {C.KOMPONENTEN[n] * t1 / Cap[i]:.2f} K, "
                                            f"mit Abgabe (1. Ordnung, Luft fest) {dT_erw:.2f} K",
              "laeufe/rom_transient/transient.csv", "CFD zwischen beiden Werten (+-10 %)",
              "OK" if 0.9 * dT_erw <= dT_cfd <= 1.1 * C.KOMPONENTEN[n] * t1 / Cap[i] else "FRAGLICH")
        # 63-%-Zeit des Bauteils CFD vs ROM, eigene Zeitkonstante C/G
        def t63(t, T):
            T = np.asarray(T)
            ziel = T[0] + 0.632 * (T[-1] - T[0])
            return float(np.interp(ziel, T, t))
        t_cfd = t63(rt["t_s"], rt[f"{n}_cfd_C"])
        t_rom = t63(rt["t_s"], rt[f"{n}_rom_C"])
        zeile(f"{n} 63-%-Zeit", f"CFD {t_cfd / 60:.0f} min, ROM {t_rom / 60:.0f} min; tau_Bauteil = C/G = {tau / 60:.0f} min "
                               f"(G {Gi:.2f} W/K)", "ergebnisse/rom/rom_transient.csv, rom_parameter.json",
              "ROM/CFD +-20 %; tau_Bauteil <= 63-%-Zeit (Luft und Wand heizen langsamer mit)",
              "OK" if abs(t_rom / t_cfd - 1) < 0.2 and tau <= t_cfd * 1.1 else "FRAGLICH")
        if n == "Schuetz" and abs(t_rom / t_cfd - 1) >= 0.2:
            text("     (Schuetz: kleine Kapazitaet, schneller Anstieg; Abweichung siehe ROM-Fehler 1,6 K)")
    tau_luft = Cap[4] / (G[4, :4].sum() + G[4, 5])
    zeile("Luft tau = C/G", f"{tau_luft / 60:.1f} min (C {Cap[4]:.0f} J/K = 1,1 kg/m3 * 1006 * V)", "rom_parameter.json",
          "Minuten (kleine Masse) -> steifes System, ROM-Eigenwert 0,2 min", "OK" if tau_luft < 600 else "FRAGLICH")
    # Steady State
    ref = pd.read_csv(HIER / "rom_export" / "referenz_python.csv")
    for z in par["stationaer"]["wert"]:
        r = ref[(ref.fall == z["name"]) & (ref.art == "stationaer")].iloc[0]
        T_rom = np.array([r[f"T_{k}_C"] for k in ["Umrichter", "Netzteil", "Schuetz", "Trafo", "Luft", "Wand"]])
        err = T_rom - (np.array(z["T_cfd_K"]) - 273.15)
        grenze = 1.5 if z["rolle"] == "Validierung" else 2.5
        zeile(f"ROM stationaer {z['name']} ({z['rolle']})", f"max |ROM-CFD| {np.max(np.abs(err)):.2f} K",
              "rom_export/referenz_python.csv vs rom_parameter.json (T_cfd aus laeufe/*/ergebnis.json)",
              f"<= {grenze} K", "OK" if np.max(np.abs(err)) <= grenze else "FRAGLICH")
    # ROM (Netz grob angepasst) gegen CFD Netz mittel
    em = lade("basis_mittel")
    r = ref[ref.fall == "basis_grob"].iloc[0]
    d = [r[f"T_{n}_C"] - em["Tmittel_C"][n] for n in BT]
    zeile("ROM (grob) gegen CFD basis_mittel", "Tmittel-Abw. " + ", ".join(f"{n} {x:+.2f}" for n, x in zip(BT, d)) + " K",
          "referenz_python.csv / laeufe/basis_mittel/ergebnis.json", "<= ROM-Fehler + Netzspannweite (~1,3 K)",
          "OK" if max(abs(x) for x in d) <= 1.3 else "FRAGLICH")
    text()


# =============================================================================================
# 5 Luefter
# =============================================================================================
def teil5(fl):
    text("5) LUEFTER-SZENARIO (35 degC) - Energiebilanz ueber die Luft")
    e = lade("sz_luefter")
    q_waende = sum(e["Q_seite_W"].values())
    rest = e["Q_quellen_W"] - q_waende
    zeile("Bilanzrest P - Q_Waende", f"{e['Q_quellen_W']:.0f} - {q_waende:.1f} = {rest:.1f} W", "laeufe/sz_luefter/ergebnis.json",
          "= Waermeabfuhr durch die Luft (Gate 6 nennt 137 W)", "OK" if abs(rest - 137) < 1.5 else "FRAGLICH",
          "Dies ist ein Rest, keine unabhaengige Messung der Luftwaerme")
    if not fl or "sz_luefter" not in fl:
        zeile("mdot cp dT", "nicht verfuegbar", "logs/audit_fluent.json fehlt", "-", "FRAGLICH")
        text()
        return
    L = fl["sz_luefter"]
    m_ein, m_aus = L["m-ein"], L["m-aus"]
    cp = L["cp_luft"]
    T_ein, T_aus = L["tm-ein"], L["tm-aus"]
    zeile("Massenerhaltung", f"m_ein {m_ein:.5f} kg/s, m_aus {m_aus:.5f} kg/s, Summe {m_ein + m_aus:+.2e}",
          "logs/audit_fluent.json (aus sz_luefter.dat.h5)", "Summe ~ 0 (< 0,5 %)",
          "OK" if abs(m_ein + m_aus) < 0.005 * abs(m_ein) else "FRAGLICH")
    zeile("Eintrittstemperatur", f"{T_ein - 273.15:.2f} degC", "logs/audit_fluent.json", "35,00 degC (Umgebung)",
          "OK" if abs(T_ein - 308.15) < 0.05 else "FRAGLICH")
    Q_luft = abs(m_ein) * cp * (T_aus - T_ein)
    zeile("P_Luft = mdot cp dT", f"{abs(m_ein):.5f} kg/s * {cp:.1f} J/(kg K) * {T_aus - T_ein:.2f} K = {Q_luft:.1f} W",
          "logs/audit_fluent.json (massengemittelte Austrittstemperatur)", f"= Bilanzrest {rest:.1f} W (+-3 %)",
          "OK" if abs(Q_luft / rest - 1) < 0.03 else "FRAGLICH",
          f"Fluent-Enthalpiefluss Oeffnungen (without-sources): {L['q-ein(without-sources)'] + L['q-aus(without-sources)']:+.1f} W")
    zeile("dT Ein/Aus", f"{T_aus - T_ein:.2f} K", "logs/audit_fluent.json",
          "Gate 6 nennt 6,6 K (dort mdot 0,0205 kg/s mit rho 1,1 statt 1,1457 gerechnet)",
          "OK" if abs((T_aus - T_ein) - 6.6) < 0.15 else "FRAGLICH",
          "Abweichung nur in der Handabschaetzung von Gate 6, nicht im CFD-Ergebnis")
    rho = 101325 / (287.05 * 308.15)
    zeile("Volumenstrom aus mdot", f"{3600 * abs(m_ein) / rho:.1f} m3/h (rho {rho:.4f} = ideales Gas 35 degC)", "logs/audit_fluent.json",
          "67,0 m3/h (ergebnis.json)", "OK" if abs(3600 * abs(m_ein) / rho - e["luefter"]["V_cfd_m3h"]) < 0.5 else "FRAGLICH")
    text()


def main():
    fl_p = HIER / "logs" / "audit_fluent.json"
    fl = json.loads(fl_p.read_text(encoding="utf-8")) if fl_p.exists() else None
    text("PLAUSIBILITAETS-AUDIT (plausibilitaet_audit.py) - prueft nur, aendert keine Ergebnisse")
    text("Format: [Urteil] Groesse: Wert | Quelle | Erwartung | Notiz.  Keine Messdaten - alle Pruefungen sind Modell gegen Modell.")
    text()
    for t in (lambda: teil1(fl), teil2, teil3, teil4, lambda: teil5(fl)):
        t()
    text("DIE 5 GROESSTEN RISIKEN FUER EINE VEROEFFENTLICHUNG")
    for r in [
        "1. Keine Messdaten: Verluste, Bauteil-k, Emissivitaeten sind Annahmen. eps 0,6 statt 0,9 -> +10..13 K,\n"
        "   k x0,5 -> Schuetz +12 K (ergebnisse/vv/vv_tabelle.csv). Absolute Temperaturen nur als Modellstudie\n"
        "   darstellen; Formulierungen wie 'validiert' vermeiden (Abgleich = Modell gegen Modell).",
        "2. Bauteile als homogene Ersatzbloecke: Tmax (z. B. Schuetz 103,8 C) ist ein Blockwert, kein Chip-/\n"
        "   Kontakttemperatur; nicht mit realen Bauteilgrenzen vergleichen. Innerer Anstieg haengt direkt an k (Annahme).",
        "3. Netzunsicherheit: GCI formal nur fuer 1 von 6 Groessen anwendbar (oszillierend, p ausserhalb 0,5..3).\n"
        "   Keine GCI-Prozentwerte als Genauigkeit zitieren; Spannweite 0,75 K plus Monitorschwankung +-0,5 K\n"
        "   nennen. Quasi-stationaere RANS-Loesung bei Ra ~ 1e9 (Uebergangsbereich) - laminar/SST nur <= 0,3 K Unterschied.",
        "4. Zirkulaere Plausibilisierung aussen: Handrechnung und CFD nutzen dieselben Aussenkorrelationen;\n"
        "   43,49 = 43,50 C beweist Konsistenz, nicht Richtigkeit des Aussenwaermeuebergangs. Das unabhaengige\n"
        "   Innenmodell (Teil 1) hat 12-23 K breite Baender.",
        "5. EM nur analytisch und ROM nur teilweise belegt: keine Feldloesung (Endeffekte, Anschluesse,\n"
        "   Gehaeusewirbelstroeme), 1D-Schienenmodell > 500 Hz ungueltig, Oberschwingungen angenommen; ROM-Fehler\n"
        "   bis 2,2 K, Fahnenkopplung nicht identifizierbar, Simulink = Python nur Konsistenz. Zusaetzlich: kleine\n"
        "   Rechenfehler in Gate-Texten moeglich (Gate 6: dT 6,6 statt 6,32 K) -> Zahlen vor Veroeffentlichung\n"
        "   nur aus den CSV/JSON uebernehmen.",
    ]:
        text(r)
    OUT.insert(2, f"Zusammenfassung: OK {ZAEHLER['OK']}, FRAGLICH {ZAEHLER['FRAGLICH']}, FEHLER {ZAEHLER['FEHLER']}")
    (HIER / "logs" / "plausibilitaet_audit.txt").write_text("\n".join(OUT) + "\n", encoding="utf-8")
    print("\n".join(OUT))


if __name__ == "__main__":
    main()
