"""
config.py - zentrale Parameter fuer das Schaltschrank-Projekt (SI-Einheiten).
Alle Werte sind ANNAHMEN fuer den Basisfall und muessen in docs/annahmen.md
begruendet werden. Aendern nur hier.

Koordinatensystem: Ursprung = Ecke hinten-links-unten des Schranks
  x = Breite (links -> rechts), y = Hoehe (unten -> oben, Gravitation -y),
  z = Tiefe (Rueckwand z=0 -> Tuer z=T)
"""

# --- Schrankgeometrie (aussen) in m ---
B = 0.80      # Breite  (x)
H = 1.20      # Hoehe   (y)
T = 0.40      # Tiefe   (z)
WAND = 0.002  # Blechdicke (m)

# --- Verlustleistung (W), Summe ueber Komponenten ---
KOMPONENTEN = {
    "Umrichter":   120.0,
    "Netzteil":     40.0,
    "Schuetz":      25.0,
    "Trafo":        65.0,
}

# --- Komponenten als Quader: Ursprungsecke (kleinstes x, y, z) und Groesse (x, y, z) in m ---
WAND_ABSTAND = 0.02   # Abstand der Komponenten zur Rueckwand (z)

BLOCK_ORIGIN = {
    "Umrichter": (0.10, 0.55, WAND_ABSTAND),
    "Netzteil":  (0.40, 0.70, WAND_ABSTAND),
    "Schuetz":   (0.40, 0.40, WAND_ABSTAND),
    "Trafo":     (0.55, 0.20, WAND_ABSTAND),
}
BLOCK_SIZE = {
    "Umrichter": (0.20, 0.40, 0.20),
    "Netzteil":  (0.15, 0.12, 0.12),
    "Schuetz":   (0.10, 0.08, 0.08),
    "Trafo":     (0.20, 0.20, 0.15),
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

# =====================================================================
# CFD (Fluent) - alles ANNAHMEN, begruendet in docs/annahmen.md
# =====================================================================
# Ersatzmaterial je Bauteil (homogener Block, isotrope effektive Waermeleitfaehigkeit)
# k in W/(m K), rho in kg/m^3, cp in J/(kg K)
BAUTEIL_MATERIAL = {
    "Umrichter": dict(k=20.0, rho=1500.0, cp=900.0),   # Alu-Kuehlkoerper + Leistungselektronik
    "Netzteil":  dict(k=5.0,  rho=1200.0, cp=900.0),   # Gehaeuse, Platine, Bauteile
    "Schuetz":   dict(k=1.5,  rho=1800.0, cp=1000.0),  # Kunststoffgehaeuse mit Cu-Kontakten
    "Trafo":     dict(k=10.0, rho=6000.0, cp=480.0),   # Blechpaket + Cu-Wicklung (Mittelwert)
}

# Schrankwand (Stahlblech, als duenne Schale mit Querleitung)
WAND_K = 45.0          # W/(m K), unlegierter Stahl
WAND_RHO = 7850.0
WAND_CP = 460.0

# Innere Emissivitaeten (S2S)
EPS_WAND_INNEN = 0.9   # lackiertes Blech innen
EPS_BAUTEIL = 0.9      # lackierte/beschichtete Bauteiloberflaechen

# Aussenseite: Boden steht auf dem Fussboden -> adiabat (wie handrechnung.py)
BODEN_ADIABAT = True
# Aussen-h wird aus den Korrelationen in handrechnung.py bestimmt (Churchill-Chu / McAdams)
# und nach dem CFD-Lauf mit den CFD-Wandtemperaturen nachgefuehrt (H_ITER_MAX Schleifen).
H_ITER_MAX = 3
H_ITER_TOL = 0.05      # relative Aenderung von h, ab der nicht mehr nachgefuehrt wird

# Turbulenzmodell Basisfall: "k-omega-sst" oder "laminar"
TURBULENZ = "k-omega-sst"
DICHTE = "incompressible-ideal-gas"   # oder "boussinesq"

# Numerik
NPROC_SOLVER = 12   # Lizenz: 16 Kerne -> "license not available", 12 ok (getestet)
ITER_MAX = 3000

# =====================================================================
# Phase 6: Szenarien (ANNAHMEN)
# =====================================================================
# Layout B: waermeempfindliche Teile unten (kuehle Luft), Verlustquellen oben -> keine Fahne trifft
# ein anderes Bauteil von unten. Abmessungen wie BLOCK_SIZE.
LAYOUTS = {
    "A": BLOCK_ORIGIN,
    "B": {
        "Umrichter": (0.55, 0.55, WAND_ABSTAND),
        "Netzteil":  (0.10, 0.20, WAND_ABSTAND),
        "Schuetz":   (0.38, 0.20, WAND_ABSTAND),
        "Trafo":     (0.10, 0.70, WAND_ABSTAND),
    },
}
# Erzwungene Kuehlung: Filterluefter unten in der Tuer (blasend), Austrittsfilter oben in der Tuer
LUEFTER_OEFFNUNG = 0.12                  # quadratische Oeffnung, Kantenlaenge m
LUEFTER_MITTE = (0.40, 0.15)             # (x, y) Mitte Eintritt in der Tuer
AUSTRITT_MITTE = (0.40, 1.05)            # (x, y) Mitte Austritt in der Tuer
# Luefterkennlinie linear: dp = DP0 * (1 - V/V0)  (typischer 120-mm-Filterluefter, Annahme)
LUEFTER_DP0 = 60.0                       # Pa bei V = 0
LUEFTER_V0 = 100.0 / 3600.0              # m^3/s freiblasend (100 m^3/h)
# Filtermatten: dp = zeta * rho/2 * v^2 (v in der Oeffnung), je Filter
FILTER_ZETA = 9.0                        # -> ca. 20 Pa bei 100 m^3/h (Annahme)

# =====================================================================
# Phase 4: EM-thermische Kopplung (ANNAHMEN, docs/annahmen.md)
# =====================================================================
# Kupfer (weichgegluehtes Cu): spez. Widerstand bei 20 degC und Temperaturbeiwert
RHO_CU_20 = 1.724e-8      # Ohm m
ALPHA_CU = 0.00393        # 1/K
CU_K, CU_RHO, CU_CP = 390.0, 8900.0, 385.0
MU0 = 4e-7 * 3.141592653589793

# --- Stromschienen L1-L3: flach, Breitseiten zueinander (Normale z), entlang x
SCHIENE_B = 0.030                 # Breite (y) m
SCHIENE_D = 0.005                 # Dicke (z) m
SCHIENE_X = (0.08, 0.72)          # Laenge 0,64 m
SCHIENE_Y0 = 1.05                 # Unterkante
SCHIENE_Z0 = {"L1": 0.050, "L2": 0.080, "L3": 0.110}   # Rueckseite der Schiene, Abstand 25 mm
I_SCHIENE = 250.0                 # A eff, Grundschwingung 50 Hz
F_NETZ = 50.0
# Oberschwingungen (Anteil der Grundschwingung), Annahme: 6-puls. Gleichrichterlast (Datacenter-Netzteile)
# h = 6k-1: Gegensystem, h = 6k+1: Mitsystem; keine durch 3 teilbaren (kein Neutralleiterstrom)
HARMONISCHE_NETZ = {1: 1.0, 5: 0.20, 7: 0.14, 11: 0.09, 13: 0.07}

# --- Trafo: Mittelfrequenz-Transformator als Fuehrungsbeispiel (Folienwicklung)
# Huellquader wie Trafo oben; in x aufgeteilt: Schenkel | Wicklung | Mittelschenkel | Wicklung | Schenkel
TRAFO_TEILE_X = [("Kern_L", 0.03), ("Wicklung_L", 0.04), ("Kern_M", 0.06), ("Wicklung_R", 0.04), ("Kern_R", 0.03)]
F_TRAFO = 20e3                    # Hz Schaltfrequenz (Annahme, DC/DC-Stufe)
FOLIE_D = 0.3e-3                  # Foliendicke m
FOLIE_B = 0.12                    # Folienbreite m
LAGEN = 8                         # Lagen (= Windungen) je Wicklung, Primaer und Sekundaer getrennt
L_WINDUNG = 0.58                  # mittlere Windungslaenge m
I_TRAFO = 60.0                    # A eff Grundschwingung (Primaer = Sekundaer, Uebersetzung 1:1)
# Rechteckstrom-Naeherung: ungerade Harmonische mit 1/n
HARMONISCHE_TRAFO = {1: 1.0, 3: 1 / 3, 5: 1 / 5, 7: 1 / 7, 9: 1 / 9, 11: 1 / 11, 13: 1 / 13}
P_KERN = 15.0                     # W Kernverluste (temperaturunabhaengig angenommen)
KERN_MAT = dict(k=20.0, rho=7650.0, cp=460.0)       # Ferrit/Blechpaket-Ersatz (isotrop)
WICKLUNG_MAT = dict(k=5.0, rho=6000.0, cp=420.0)    # Cu-Folie + Isolation, quer zur Folie schlecht leitend
SCHIENE_MAT = dict(k=CU_K, rho=CU_RHO, cp=CU_CP)
EM_TOL_P = 0.005                  # Abbruch Kopplung: rel. Aenderung der Verluste
EM_TOL_T = 0.1                    # K
EM_ITER_MAX = 8
