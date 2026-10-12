"""
handrechnung.py - Phase 1: physikalische Abschaetzung des Schrank-Temperaturanstiegs.

Verfahren: Waermebilanz der Schrankoberflaeche mit
  - natuerlicher Konvektion an senkrechten Waenden (Churchill-Chu) bzw. Deckel,
  - Strahlung an die Umgebung,
  - gleichmaessiger Wandtemperatur (grobe Naeherung, ignoriert Schichtung innen).
Ergebnis ist eine ABSCHAETZUNG zur Plausibilisierung der CFD, kein Normnachweis.

IEC 60890 / DIN EN 61439-1: siehe Abschnitt unten. Normkoeffizienten sind
bewusst NICHT eingetragen (Originaltext erforderlich).
"""

import math
from config import *

SIGMA = 5.670374e-8


def kelvin(c):
    return c + 273.15


def oberflaechen():
    """Frei abstrahlende/abkuehlende Flaechen eines freistehenden Schranks (m^2).
    Boden wird als nicht beteiligt angenommen."""
    seiten = 2 * H * T            # 2 Seitenwaende
    front_rueck = 2 * H * B       # Front + Rueck
    deckel = B * T
    return {"senkrecht": seiten + front_rueck, "deckel": deckel}


def nusselt_senkrecht(Ra, Pr):
    """Churchill-Chu fuer senkrechte Platte (ueber gesamten Ra-Bereich)."""
    return (0.825 + 0.387 * Ra ** (1 / 6) / (1 + (0.492 / Pr) ** (9 / 16)) ** (8 / 27)) ** 2


def nusselt_deckel_heiss_nach_oben(Ra):
    """Waagerechte Platte, warme Seite oben (McAdams)."""
    if Ra < 1e7:
        return 0.54 * Ra ** 0.25
    return 0.15 * Ra ** (1 / 3)


def h_konvektion(T_wand_c, lage, T_umg_c=None):
    """Aeusserer konvektiver Waermeuebergangskoeffizient W/(m^2 K) fuer eine Schrankflaeche.
    lage: 'senkrecht' (Churchill-Chu, L = H) oder 'deckel' (McAdams, warme Seite oben, L = A/U).
    Stoffwerte konstant aus config.py (ca. 45 degC). Fuer dT <= 0.1 K wird dT = 0.1 K gesetzt."""
    T_umg_c = T_UMGEBUNG if T_umg_c is None else T_umg_c
    Tw, Ta = kelvin(T_wand_c), kelvin(T_umg_c)
    dT = max(Tw - Ta, 0.1)
    beta = 1.0 / ((Tw + Ta) / 2)
    if lage == "senkrecht":
        L = H
        Ra = G * beta * dT * L ** 3 / (LUFT_NU * LUFT_ALPHA)
        return nusselt_senkrecht(Ra, LUFT_PR) * LUFT_LAMBDA / L
    if lage == "deckel":
        L = (B * T) / (2 * (B + T))
        Ra = G * beta * dT * L ** 3 / (LUFT_NU * LUFT_ALPHA)
        return nusselt_deckel_heiss_nach_oben(Ra) * LUFT_LAMBDA / L
    raise ValueError(lage)


def waermeabgabe(T_wand_c, T_umg_c=None, eps=None, hfaktor=1.0):
    """Abgegebene Leistung (W) bei gegebener Wandtemperatur (degC)."""
    T_umg_c = T_UMGEBUNG if T_umg_c is None else T_umg_c
    eps = EMISSIVITAET if eps is None else eps
    A = oberflaechen()
    Tw, Ta = kelvin(T_wand_c), kelvin(T_umg_c)
    dT = Tw - Ta
    if dT <= 0:
        return 0.0, {}

    beta = 1.0 / ((Tw + Ta) / 2)

    # senkrechte Waende: charakteristische Laenge = Hoehe
    Ra_v = G * beta * dT * H ** 3 / (LUFT_NU * LUFT_ALPHA)
    Nu_v = nusselt_senkrecht(Ra_v, LUFT_PR)
    h_v = Nu_v * LUFT_LAMBDA / H

    # Deckel: charakteristische Laenge = A/U
    L_d = (B * T) / (2 * (B + T))
    Ra_d = G * beta * dT * L_d ** 3 / (LUFT_NU * LUFT_ALPHA)
    Nu_d = nusselt_deckel_heiss_nach_oben(Ra_d)
    h_d = Nu_d * LUFT_LAMBDA / L_d

    A_ges = A["senkrecht"] + A["deckel"]
    h_v, h_d = hfaktor * h_v, hfaktor * h_d
    q_konv = h_v * A["senkrecht"] * dT + h_d * A["deckel"] * dT
    q_str = eps * SIGMA * (Tw ** 4 - Ta ** 4) * A_ges

    info = {
        "Ra_senkrecht": Ra_v, "Nu_senkrecht": Nu_v, "h_konv_senkrecht": h_v,
        "Ra_deckel": Ra_d, "h_konv_deckel": h_d,
        "Q_konvektion": q_konv, "Q_strahlung": q_str,
        "A_gesamt": A_ges,
    }
    return q_konv + q_str, info


def loese_wandtemperatur(Q, T_umg_c=None, eps=None, hfaktor=1.0, tol=1e-6):
    """Bisektion: Wandtemperatur, bei der die abgegebene Leistung gleich Q ist."""
    T_umg_c = T_UMGEBUNG if T_umg_c is None else T_umg_c
    lo, hi = T_umg_c, T_umg_c + 200.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        q, _ = waermeabgabe(mid, T_umg_c, eps, hfaktor)
        if q < Q:
            lo = mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi)


def main():
    Q = sum(KOMPONENTEN.values())
    Tw = loese_wandtemperatur(Q)
    q, info = waermeabgabe(Tw)

    print("=" * 60)
    print("Handrechnung Schaltschrank (Abschaetzung)")
    print("=" * 60)
    print(f"Abmessungen B x H x T : {B} x {H} x {T} m")
    print(f"Verlustleistung gesamt: {Q:.1f} W")
    for k, v in KOMPONENTEN.items():
        print(f"   {k:12s}: {v:6.1f} W")
    print(f"Umgebung              : {T_UMGEBUNG:.1f} degC")
    print(f"Flaeche gesamt (ohne Boden): {info['A_gesamt']:.3f} m^2")
    print()
    print(f"Ra senkrecht          : {info['Ra_senkrecht']:.3e}")
    print(f"Nu senkrecht          : {info['Nu_senkrecht']:.1f}")
    print(f"h konv. senkrecht     : {info['h_konv_senkrecht']:.2f} W/(m^2 K)")
    print(f"h konv. Deckel        : {info['h_konv_deckel']:.2f} W/(m^2 K)")
    reg = "turbulent (Uebergangsbereich/turbulent)" if info["Ra_senkrecht"] > 1e9 else "laminar"
    print(f"Stroemung aussen      : {reg}")
    print()
    print(f"Wandtemperatur (gleichmaessig angenommen): {Tw:.1f} degC "
          f"(Anstieg {Tw - T_UMGEBUNG:.1f} K)")
    print(f"   davon Konvektion   : {info['Q_konvektion']:.1f} W")
    print(f"   davon Strahlung    : {info['Q_strahlung']:.1f} W "
          f"({100 * info['Q_strahlung'] / Q:.0f} % der Gesamtleistung)")
    print()
    print("Hinweis: Die Innenlufttemperatur liegt im Betrieb DEUTLICH ueber der")
    print("Wandtemperatur (Schichtung, Waermeuebergang innen, Blechwiderstand).")
    print("Der Innenanstieg wird erst durch die CFD aufgeloest; dieser Wert ist")
    print("eine Untergrenze fuer die mittlere Innenluft.")

    # ---------------------------------------------------------------
    # IEC 60890 - Platzhalter
    # ---------------------------------------------------------------
    print()
    print("-" * 60)
    print("IEC 60890 / DIN EN 61439-1 : noch nicht berechnet.")
    print("Die Normkoeffizienten (Oberflaechenfaktoren, Einbaufaktoren, Verteilungs-")
    print("Faktoren, Exponenten) muessen aus dem Originaltext uebernommen werden.")
    print("Bis dahin dient obige physikalische Abschaetzung als Plausibilitaetscheck.")
    print("-" * 60)


if __name__ == "__main__":
    main()
