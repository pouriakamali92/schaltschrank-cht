"""
config.py - zentrale Parameter fuer das Schaltschrank-Projekt (SI-Einheiten).
Alle Werte sind ANNAHMEN fuer den Basisfall und muessen in docs/annahmen.md
begruendet werden. Aendern nur hier.
"""

# --- Schrankgeometrie (aussen) in m ---
B = 0.80      # Breite
H = 1.20      # Hoehe
T = 0.40      # Tiefe
WAND = 0.002  # Blechdicke (m)

# --- Verlustleistung (W), Summe ueber Komponenten ---
KOMPONENTEN = {
    "Umrichter":   120.0,
    "Netzteil":     40.0,
    "Schuetz":      25.0,
    "Trafo":        65.0,
}

# --- Umgebung ---
T_UMGEBUNG = 35.0   # degC
T_MAX_INNEN = 55.0  # degC, typ. Auslegungsgrenze Innenluft (Annahme)

# --- Oberflaeche ---
EMISSIVITAET = 0.9  # lackiertes Blech (Annahme)

# --- Stoffwerte Luft (bei ca. 45 degC, 1 atm; Naeherungswerte) ---
LUFT_NU = 1.75e-5      # kinematische Viskositaet m^2/s
LUFT_ALPHA = 2.45e-5   # Temperaturleitfaehigkeit m^2/s
LUFT_LAMBDA = 0.0275   # Waermeleitfaehigkeit W/(m K)
LUFT_PR = 0.71

G = 9.81
