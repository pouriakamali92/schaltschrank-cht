"""
em_modell.py - reduziertes analytisches EM-Modell der ohmschen Verluste (Phase 4).

KEINE Feldloesung mit FEM (Ansys Maxwell ist auf diesem Rechner nicht installiert).
Modelle (alle Formeln hier im Code nachvollziehbar):
  rho(T)   = rho20 * (1 + alpha * (T - 20 degC))
  delta    = sqrt(rho / (pi f mu0))                      Skineindringtiefe
  1D-Platte: magnetoquasistatische Diffusion in einer unendlich breiten Platte der Dicke t mit
           komplexen Oberflaechenfeldern Ha, Hb (Effektivwert-Zeiger):
             H(x) = [Ha sinh(k (t-x)) + Hb sinh(k x)] / sinh(k t),  k = (1+j)/delta,  J = dH/dx
             P/A  = rho * Integral |J|^2 dx                (A = Plattenflaeche)
           -> Skin- UND Proximity-Effekt, sobald die Felder der Nachbarleiter in Ha, Hb stecken.
  Dowell:  F_R = Delta [ s1 + 2/3 (m^2 - 1) s2 ],  Delta = sqrt(eta) h / delta,
           s1 = (sinh 2D + sin 2D)/(cosh 2D - cos 2D),  s2 = (sinh D - sin D)/(cosh D + cos D)
           (Folienwicklung mit m Lagen; eta = Folienbreite / Fensterhoehe)
  2D-Filament-Methode (Pruefung der 1D-Naeherung fuer die Schienen): Leiterquerschnitte in
           Teilleiter zerlegt, Teilinduktivitaeten je Laenge M_ij = -mu0/(2 pi) ln d_ij
           (Selbstterm: mittlerer geometrischer Abstand des Rechtecks ~ 0.2235 (a+b)),
           gleiche Spannung je Schiene, vorgegebener Schienenstrom -> Stromverteilung -> Verluste.

Aufruf: python em_modell.py   -> ergebnisse/em/*.png, *.csv, Konsolenausgabe
"""
import cmath
import csv
import math
from pathlib import Path

import numpy as np

import config as C

HIER = Path(__file__).resolve().parent
ZIEL = HIER / "ergebnisse" / "em"
PHASEN = list(C.SCHIENE_Z0)            # L1, L2, L3 (Reihenfolge entlang z)


# ---------------------------------------------------------------------------
# Grundgroessen
# ---------------------------------------------------------------------------
def rho_cu(T_c):
    return C.RHO_CU_20 * (1.0 + C.ALPHA_CU * (T_c - 20.0))


def skintiefe(f, T_c):
    return math.sqrt(rho_cu(T_c) / (math.pi * f * C.MU0))


def platte_verlust(Ha, Hb, t, f, T_c, n=4001):
    """Verlust je Plattenflaeche (W/m^2) einer Platte der Dicke t mit Oberflaechenfeldern Ha, Hb
    (A/m, komplexe Effektivwert-Zeiger). Gleichstromgrenzfall fuer f -> 0 analytisch."""
    rho = rho_cu(T_c)
    if f <= 0:
        return rho * abs(Hb - Ha) ** 2 / t
    d = skintiefe(f, T_c)
    k = (1 + 1j) / d
    x = np.linspace(0.0, t, n)
    sh = np.sinh(k * t)
    J = (-Ha * k * np.cosh(k * (t - x)) + Hb * k * np.cosh(k * x)) / sh
    return rho * np.trapezoid(np.abs(J) ** 2, x)


def dowell_FR(delta_h, m):
    """Dowell-Faktor R_ac/R_dc fuer m Lagen, normierte Lagendicke delta_h = sqrt(eta) h / delta."""
    D = delta_h
    s1 = (math.sinh(2 * D) + math.sin(2 * D)) / (math.cosh(2 * D) - math.cos(2 * D))
    s2 = (math.sinh(D) - math.sin(D)) / (math.cosh(D) + math.cos(D))
    return D * (s1 + 2.0 / 3.0 * (m * m - 1) * s2)


# ---------------------------------------------------------------------------
# Stromschienen
# ---------------------------------------------------------------------------
def schienen_stroeme(h, I1):
    """Effektivwert-Zeiger der drei Phasenstroeme der h-ten Harmonischen.
    h = 6k+1 Mitsystem, h = 6k-1 Gegensystem (symmetrische Last, keine Nullkomponente)."""
    folge = 1 if h % 6 == 1 else -1
    return [I1 * cmath.exp(-1j * folge * 2 * math.pi / 3 * i) for i in range(3)]


def schienen_1d(f, T_c, stroeme, w=C.SCHIENE_B, t=C.SCHIENE_D):
    """Verlust je Laenge (W/m) jeder Schiene, 1D-Modell (Schienen als unendlich breite Platten,
    Breitseiten zueinander). Feld zwischen den Schienen aus dem Durchflutungsgesetz:
    H links von Schiene 1 = 0, H springt an jeder Schiene um I/w."""
    P = []
    Ha = 0j
    for I in stroeme:
        Hb = Ha + I / w
        P.append(platte_verlust(Ha, Hb, t, f, T_c) * w)
        Ha = Hb
    return P


def einzelschiene_1d(f, T_c, I=1.0, w=C.SCHIENE_B, t=C.SCHIENE_D):
    """Isolierte Schiene (nur Skineffekt): symmetrische Felder -I/(2w), +I/(2w)."""
    return platte_verlust(-I / (2 * w), I / (2 * w), t, f, T_c) * w


def schienen_2d(f, T_c, stroeme, ny=36, nz=9):
    """2D-Filament-Methode: Verlust je Laenge (W/m) jeder Schiene."""
    w, t = C.SCHIENE_B, C.SCHIENE_D
    a, b = w / ny, t / nz
    pos, bar = [], []
    for k, ph in enumerate(PHASEN):
        z0 = C.SCHIENE_Z0[ph]
        for iy in range(ny):
            for iz in range(nz):
                pos.append(((iy + 0.5) * a, z0 + (iz + 0.5) * b))
                bar.append(k)
    pos = np.array(pos)
    N = len(pos)
    d = np.hypot(pos[:, None, 0] - pos[None, :, 0], pos[:, None, 1] - pos[None, :, 1])
    np.fill_diagonal(d, 0.2235 * (a + b))
    M = -C.MU0 / (2 * math.pi) * np.log(d)
    R = rho_cu(T_c) / (a * b)
    w_om = 2 * math.pi * f
    # Unbekannte: N Teilstroeme + 3 Spannungen je Laenge
    A = np.zeros((N + 3, N + 3), dtype=complex)
    rhs = np.zeros(N + 3, dtype=complex)
    A[:N, :N] = 1j * w_om * M + R * np.eye(N)
    for i in range(N):
        A[i, N + bar[i]] = -1.0
    for k in range(3):
        A[N + k, :N] = [1.0 if bar[i] == k else 0.0 for i in range(N)]
        rhs[N + k] = stroeme[k]
    x = np.linalg.solve(A, rhs)
    i_f = x[:N]
    return [float(R * np.sum(np.abs(i_f[np.array(bar) == k]) ** 2)) for k in range(3)]


def schienen_verluste(T_schiene, modell="1d"):
    """Verluste (W) je Schiene inkl. Oberschwingungen bei Schienentemperaturen T_schiene (dict)."""
    L = C.SCHIENE_X[1] - C.SCHIENE_X[0]
    out = {ph: 0.0 for ph in PHASEN}
    for h, anteil in C.HARMONISCHE_NETZ.items():
        f = h * C.F_NETZ
        I = schienen_stroeme(h, C.I_SCHIENE * anteil)
        for k, ph in enumerate(PHASEN):
            # Temperatur je Schiene: Modell mit der Temperatur der jeweiligen Schiene auswerten
            P = (schienen_1d if modell == "1d" else schienen_2d)(f, T_schiene[ph], I)
            out[ph] += P[k] * L
    return out


def schienen_rdc(T_c):
    L = C.SCHIENE_X[1] - C.SCHIENE_X[0]
    return rho_cu(T_c) * L / (C.SCHIENE_B * C.SCHIENE_D)


# ---------------------------------------------------------------------------
# Trafo-Wicklung (Mittelfrequenz, Folie)
# ---------------------------------------------------------------------------
FENSTER_H = C.BLOCK_SIZE["Trafo"][1]
ETA = C.FOLIE_B / FENSTER_H


def wicklung_rdc(T_c):
    """Gleichstromwiderstand Primaer + Sekundaer (je LAGEN Windungen)."""
    return 2 * rho_cu(T_c) * C.LAGEN * C.L_WINDUNG / (C.FOLIE_D * C.FOLIE_B)


def wicklung_FR(f, T_c, m=C.LAGEN):
    return dowell_FR(math.sqrt(ETA) * C.FOLIE_D / skintiefe(f, T_c), m)


def wicklung_verlust(T_c):
    """Cu-Verluste (W) der Wicklung inkl. Harmonischer des Rechteckstroms."""
    Rdc = wicklung_rdc(T_c)
    return sum(Rdc * (C.I_TRAFO * a) ** 2 * wicklung_FR(h * C.F_TRAFO, T_c)
               for h, a in C.HARMONISCHE_TRAFO.items())


def dowell_pruefung(f=C.F_TRAFO, T_c=80.0, m=C.LAGEN):
    """Dowell-Formel gegen die schichtweise 1D-Plattenloesung (Lage l: Ha=(l-1)I/b, Hb=l I/b)."""
    b = C.FOLIE_B / math.sqrt(ETA)     # effektive Feldbreite (Porositaet)
    h_eff = C.FOLIE_D
    # Mit Porositaet: Platte der Dicke h, Leitfaehigkeit*eta -> Delta = sqrt(eta) h/delta.
    # Direkter Vergleich ohne Porositaet (eta = 1):
    d = skintiefe(f, T_c)
    P = sum(platte_verlust((l - 1) * 1.0, l * 1.0, h_eff, f, T_c) for l in range(1, m + 1))
    P_dc = m * rho_cu(T_c) * 1.0 / h_eff
    return P / P_dc, dowell_FR(h_eff / d, m)


# ---------------------------------------------------------------------------
# Verlustverteilung auf die CFD-Zonen
# ---------------------------------------------------------------------------
def zonen_verluste(T):
    """T: dict mit 'L1','L2','L3','wicklung' (degC, volumengemittelt aus der CFD).
    Rueckgabe: Verlust (W) je Fluent-Zone."""
    Ps = schienen_verluste({ph: T[ph] for ph in PHASEN}, "2d")   # 2D-Filament = Referenz
    P_cu = wicklung_verlust(T["wicklung"])
    out = {f"schiene_{ph.lower()}-schiene_{ph.lower()}": Ps[ph] for ph in PHASEN}
    out["trafo_wicklung_l-trafo_wicklung_l"] = P_cu / 2
    out["trafo_wicklung_r-trafo_wicklung_r"] = P_cu / 2
    kern = {n: b for n, b in C.TRAFO_TEILE_X if n.startswith("Kern")}
    for n, b in kern.items():
        out[f"trafo_{n.lower()}-trafo_{n.lower()}"] = C.P_KERN * b / sum(kern.values())
    return out


# ---------------------------------------------------------------------------
# Auswertung / Bilder
# ---------------------------------------------------------------------------
def main():
    import plotstil as PS
    PS.stil()
    plt = PS.plt
    ZIEL.mkdir(parents=True, exist_ok=True)
    T_ref = 80.0
    print(f"Skintiefe Cu bei {T_ref:.0f} degC: 50 Hz {1e3 * skintiefe(50, T_ref):.2f} mm, "
          f"20 kHz {1e3 * skintiefe(20e3, T_ref):.3f} mm")
    r1, r2 = dowell_pruefung()
    print(f"Pruefung Dowell (m={C.LAGEN}, 20 kHz, eta=1): Plattenloesung {r1:.4f}, Dowell-Formel {r2:.4f}")

    # --- Schienen R_ac/R_dc ueber f (Grundschwingung symmetrisch, Mitsystem)
    fs = np.logspace(0, 4, 41)
    zeilen = []
    for f in fs:
        I = schienen_stroeme(1, 1.0)
        p1 = schienen_1d(f, T_ref, I)
        p2 = schienen_2d(f, T_ref, I, 48, 12)
        pe = einzelschiene_1d(f, T_ref)
        pdc = rho_cu(T_ref) / (C.SCHIENE_B * C.SCHIENE_D)
        zeilen.append([f, pe / pdc] + [x / pdc for x in p1] + [x / pdc for x in p2])
    z = np.array(zeilen)
    with open(ZIEL / "schienen_rac_rdc.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["f_Hz", "einzeln_1D", "L1_1D", "L2_1D", "L3_1D", "L1_2D", "L2_2D", "L3_2D"])
        wr.writerows(np.round(z, 6).tolist())
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.semilogx(z[:, 0], z[:, 1], color=PS.NEUTRAL, label="einzelne Schiene (nur Skin), 1D")
    ax.semilogx(z[:, 0], z[:, 2], color=PS.KATEGORIE[0], label="aeussere Schiene L1, 1D")
    ax.semilogx(z[:, 0], z[:, 3], color=PS.KATEGORIE[1], label="mittlere Schiene L2, 1D")
    ax.semilogx(z[:, 0], z[:, 5], color=PS.KATEGORIE[0], ls="--", label="L1, 2D-Filament")
    ax.semilogx(z[:, 0], z[:, 6], color=PS.KATEGORIE[1], ls="--", label="L2, 2D-Filament")
    for h in C.HARMONISCHE_NETZ:
        ax.axvline(h * C.F_NETZ, color=PS.GITTER, lw=1, zorder=0)
    ax.text(C.F_NETZ * 0.9, 9.0, "50 Hz und\nHarmonische\n(5., 7., 11., 13.)", fontsize=8,
            color=PS.TEXT2, ha="right")
    ax.set_xlabel("Frequenz in Hz")
    ax.set_ylabel("R_ac / R_dc")
    ax.set_title(f"Cu-Stromschienen {1e3 * C.SCHIENE_B:.0f} x {1e3 * C.SCHIENE_D:.0f} mm, Abstand "
                 f"{1e3 * (C.SCHIENE_Z0['L2'] - C.SCHIENE_Z0['L1'] - C.SCHIENE_D):.0f} mm, 3-phasig, "
                 f"{T_ref:.0f} °C\nSkin- und Proximity-Effekt (analytisch 1D und 2D-Filament)", fontsize=9)
    ax.legend(fontsize=8)
    fig.savefig(ZIEL / "schienen_rac_rdc.png")
    plt.close(fig)

    # --- Trafo: Dowell ueber f fuer mehrere Lagenzahlen
    fs2 = np.logspace(2, 6, 61)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    zeilen = []
    for i, m in enumerate((1, 2, 4, 8)):
        F = [wicklung_FR(f, T_ref, m) for f in fs2]
        zeilen.append(F)
        ax.loglog(fs2, F, color=PS.KATEGORIE[i], label=f"m = {m} Lagen" + (" (Auslegung)" if m == C.LAGEN else ""))
    ax.axvline(C.F_TRAFO, color=PS.TEXT2, lw=1, ls=":")
    ax.text(C.F_TRAFO * 1.08, 1.05, f"{C.F_TRAFO / 1e3:.0f} kHz", fontsize=8, color=PS.TEXT2)
    ax.set_xlabel("Frequenz in Hz")
    ax.set_ylabel("R_ac / R_dc (Dowell)")
    ax.set_title(f"Folienwicklung {1e3 * C.FOLIE_D:.1f} mm Cu, Porositaet {ETA:.2f}, {T_ref:.0f} °C", fontsize=9)
    ax.legend(fontsize=8)
    fig.savefig(ZIEL / "trafo_dowell.png")
    plt.close(fig)
    with open(ZIEL / "trafo_dowell.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["f_Hz", "m1", "m2", "m4", "m8"])
        for j, f in enumerate(fs2):
            wr.writerow([f] + [round(zeilen[i][j], 6) for i in range(4)])

    # --- Betriebspunkt
    T = {"L1": T_ref, "L2": T_ref, "L3": T_ref, "wicklung": T_ref}
    p1 = schienen_verluste(T, "1d")
    p2 = schienen_verluste(T, "2d")
    Pdc = {ph: schienen_rdc(T_ref) * C.I_SCHIENE ** 2 for ph in PHASEN}
    print(f"Schienen bei {T_ref:.0f} degC (I1 = {C.I_SCHIENE:.0f} A + Harmonische):")
    for ph in PHASEN:
        print(f"  {ph}: P_dc(I1) {Pdc[ph]:.2f} W, P_ac 1D {p1[ph]:.2f} W, 2D {p2[ph]:.2f} W "
              f"(Abweichung 1D/2D {100 * (p1[ph] / p2[ph] - 1):+.1f} %)")
    Rdc = wicklung_rdc(T_ref)
    I_eff = C.I_TRAFO * math.sqrt(sum(a * a for a in C.HARMONISCHE_TRAFO.values()))
    print(f"Wicklung: R_dc {1e3 * Rdc:.2f} mOhm, I_eff {I_eff:.1f} A, P_dc {Rdc * I_eff ** 2:.1f} W, "
          f"P_ac {wicklung_verlust(T_ref):.1f} W, F_R(20 kHz) {wicklung_FR(C.F_TRAFO, T_ref):.2f}")
    for Tc in (40, 80, 120):
        print(f"  T = {Tc:3d} degC: Wicklung {wicklung_verlust(Tc):6.2f} W, R_dc-Faktor "
              f"{rho_cu(Tc) / rho_cu(20):.3f}, F_R(20 kHz) {wicklung_FR(C.F_TRAFO, Tc):.3f}")
    print("Zonenverluste bei 80 degC:", {k: round(v, 2) for k, v in zonen_verluste(T).items()})


if __name__ == "__main__":
    main()
