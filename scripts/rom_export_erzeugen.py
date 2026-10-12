"""
rom_export_erzeugen.py - erzeugt den Ordner rom_export/ fuer den Simulink-Lauf auf dem Laptop.

Inhalt:
  rom_export/rom_parameter.json   alle ROM-Parameter mit Einheit und Quelle (aus ergebnisse/rom/rom_parameter.json,
                                  config.py und laeufe/*/ergebnis.json; NICHT neu angepasst)
  rom_export/referenz_python.csv  Python-ROM (gleiche Gleichungen wie rom.py) fuer dieselben Faelle wie Simulink
Pruefungen (Ausgabe nach logs/rom_export.txt):
  - JSON neu laden, Groessen/Vorzeichen/Endlichkeit, Kapazitaeten = rho c V, Betriebspunkt reproduziert
  - Energiebilanz des ROM im stationaeren Zustand
  - Python-ROM gegen CFD (stationaer alle 6 Faelle, transient 8 h)

Aufruf: python rom_export_erzeugen.py
"""
import csv
import json
import re
import time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

import config as C
import rom as R

HIER = Path(__file__).resolve().parent
EXPORT = HIER / "rom_export"
QUELLE_PAR = HIER / "ergebnisse" / "rom" / "rom_parameter.json"
T_STAT_H = 48.0          # Simulationsdauer fuer "stationaer" (langsamste Zeitkonstante ~2,1 h)
DT_AUSGABE = 300.0       # Ausgaberaster transient = Raster des Fluent-Laufs (s)


def netz_aus_json(p):
    theta = np.r_[np.log(np.maximum(np.r_[p["gL"], p["gW"], p["gLW"]], 1e-300)), p["s"]]
    return R.Netz(theta)


def hA(h):
    return float(sum(h[s] * R.SEITE_A[s] for s in R.SEITE_A))


def wert(v, einheit, quelle):
    return {"wert": v, "einheit": einheit, "quelle": quelle}


def g_matrix(net, T, hAw, T_U, eps):
    """Sekanten-Leitwerte G_ij = q_ij / (T_i - T_j) am Betriebspunkt (W/K), symmetrisch, 0 = keine Verbindung.
    Knoten 0..5 wie KNOTEN; die Umgebung ist ein Randknoten (eigene Spalte G_umg)."""
    G = np.zeros((6, 6))
    for i in range(4):
        dL = T[i] - T[4]
        G[i, 4] = G[4, i] = net.gL[i] * (dL * dL + 1) ** (1 / 6)
        G[i, 5] = G[5, i] = net.gW[i] * R.SIGMA * (T[i] ** 2 + T[5] ** 2) * (T[i] + T[5])
    dLW = T[4] - T[5]
    G[4, 5] = G[5, 4] = net.gLW * (dLW * dLW + 1) ** (1 / 6)
    G_umg = hAw + eps * R.SIGMA * R.A_AB * (T[5] ** 2 + T_U ** 2) * (T[5] + T_U)
    return G, G_umg


# Statische Pruefung rom_simulink.m: nur Bloecke der Simulink-Basisbibliothek, keine Toolbox-Funktionen,
# alle verwendeten Funktionen spaetestens in R2022a vorhanden (Einfuehrungs-Release lt. MATLAB-Doku).
FKT_RELEASE = {"jsondecode": "R2016b", "tiledlayout": "R2019b", "nexttile": "R2019b", "exportgraphics": "R2020a",
               "datetime": "R2014b", "readtable": "R2013b", "strjoin": "R2013a", "ode15s": "vor R2006a",
               "odeset": "vor R2006a", "str2func": "vor R2006a", "interp1": "vor R2006a", "license": "vor R2006a",
               "new_system": "Simulink", "add_block": "Simulink", "add_line": "Simulink", "set_param": "Simulink",
               "sim": "Simulink", "save_system": "Simulink", "close_system": "Simulink", "bdIsLoaded": "Simulink"}
BLOECKE_BASIS = {"simulink/Sources/Constant", "simulink/Signal Routing/Mux", "simulink/User-Defined Functions/Fcn",
                 "simulink/Continuous/Integrator", "simulink/Sinks/To Workspace"}
TOOLBOX_HINWEISE = [r"\bfindop\b", r"\blinearize\b", r"\boperpoint\b", r"\bslLinearizer\b", r"\bsltest",
                    r"Simulink\.sdi", r"Stateflow", r"\bsfroot\b", r"MATLAB Function", r"Simscape", r"\bssc_",
                    r"\boptimset\b", r"\blsqnonlin\b", r"\bfsolve\b", r"\bss\s*\(", r"\btf\s*\("]


def pruefe_matlab_skript(pfad):
    src = pfad.read_text(encoding="utf-8")
    code = "\n".join(z.split("%")[0] for z in src.splitlines())          # ohne Kommentare
    bloecke = set(re.findall(r"add_block\('([^']+)'", code))
    fremd = sorted(b for b in bloecke if b not in BLOECKE_BASIS)
    toolbox = [t for t in TOOLBOX_HINWEISE if re.search(t, code)]
    nicht_ascii = sum(ord(c) > 127 for c in src)
    benutzt = sorted(f for f in FKT_RELEASE if re.search(rf"\b{re.escape(f)}\s*\(", code))
    # new_system ist ein eingebauter Simulink-Befehl: exist('new_system','file') liefert 0 auch MIT Simulink
    # (Fehler im ersten Laptop-Lauf, R2024a) -> nur exist('new_system') ist korrekt.
    exist_ok = (re.search(r"exist\('new_system'\)", code) is not None
                and re.search(r"exist\('new_system'\s*,", code) is None)
    return dict(bloecke=sorted(bloecke), fremd=fremd, toolbox=toolbox, nicht_ascii=nicht_ascii, exist_ok=exist_ok,
                funktionen=[f"{f} ({FKT_RELEASE[f]})" for f in benutzt])


def main():
    EXPORT.mkdir(exist_ok=True)
    log = []
    p = json.loads(QUELLE_PAR.read_text(encoding="utf-8"))
    net = netz_aus_json(p)
    Cap = np.array(p["C"])
    KN = p["knoten"]

    # ---------------- Faelle
    stat = []
    for fall in R.FIT + R.VAL:
        x = R.daten(fall)
        stat.append(dict(
            name=fall, rolle="Anpassung" if fall in R.FIT else "Validierung",
            P_W=x["P"].tolist(), T_umg_K=float(x["T_U"]), hA_aussen_W_per_K=hA(x["h"]),
            T0_K=[float(x["T_U"])] * 6, t_ende_s=T_STAT_H * 3600,
            T_cfd_K=x["T"].tolist(),
            quelle=f"laeufe/{fall}/ergebnis.json (Fluent stationaer, Netz grob): P, T_umg, h_aussen, "
                   "T_cfd = Volumenmittel der Bauteile, Luft-Mittel, flaechengewichtete Wand"))
    tr_json = json.loads((HIER / "laeufe" / "rom_transient" / "ergebnis.json").read_text(encoding="utf-8"))
    import pandas as pd
    d = pd.read_csv(HIER / "laeufe" / "rom_transient" / "transient.csv")
    t_end = float(d["t"].iloc[-1])
    basis = R.daten("basis_grob")
    trans = dict(
        name="rom_transient", rolle="Validierung transient",
        P_W=basis["P"].tolist(), T_umg_K=float(basis["T_U"]), hA_aussen_W_per_K=hA(tr_json["h_aussen"]),
        T0_K=[float(basis["T_U"])] * 6, t_ende_s=t_end, dt_ausgabe_s=DT_AUSGABE,
        quelle="laeufe/rom_transient (Fluent transient, dt 30 s, 8 h): Start alles auf T_umg, "
               "alle Verluste ab t = 0 eingeschaltet; Verlauf in laeufe/rom_transient/transient.csv")

    G, G_umg = g_matrix(net, np.array(p["T_ss_basis"]), hA(p["h_basis"]), p["T_U_basis"], p["eps_aussen"])
    with np.errstate(divide="ignore"):
        Rm = np.where(G > 0, 1.0 / np.where(G > 0, G, 1.0), np.nan)

    q_anp = "Least-Squares-Anpassung an 5 stationaere Fluent-Laeufe (rom.py, ergebnisse/rom/rom_parameter.json)"
    par = {
        "_info": ("Thermisches ROM (RC-Netzwerk) Schaltschrank, Selbststudienprojekt. Alle Eingaben sind Annahmen "
                  "(docs/annahmen.md), Referenz ist Fluent-CFD, keine Messdaten. Erzeugt von rom_export_erzeugen.py "
                  f"am {time.strftime('%Y-%m-%d %H:%M')}."),
        "gleichungen": {
            "Bauteil_i": "C_i dT_i/dt = P_i - qL_i - qW_i",
            "Luft": "C_L dT_L/dt = sum(qL_i) - qLW",
            "Wand": "C_W dT_W/dt = sum(qW_i) + qLW - hA (T_W - T_umg) - eps sigma A_ab (T_W^4 - T_umg^4)",
            "qL_i": "gL_i ((T_i - T_L)^2 + 1 K^2)^(1/6) (T_i - T_L)   [freie Konvektion, h ~ dT^(1/3), geglaettet]",
            "qW_i": "gW_i sigma (T_i^4 - T_W^4)   [Strahlung Bauteil -> Wand]",
            "qLW": "gLW ((T_L - T_W)^2 + 1 K^2)^(1/6) (T_L - T_W)",
            "Tmax_i": "T_i + r_int_i P_i (Nachrechnung, nicht Teil der Simulation)",
        },
        "knoten": wert(KN, "-", "rom.py: 4 Bauteil-Mitteltemperaturen, Luft-Mittel, Wand-Mittel"),
        "temperaturen": "alle Temperaturen in K, Ausgabe-CSV in degC",
        "C": wert(p["C"], "J/K", "rho cp V aus config.py (BAUTEIL_MATERIAL, BLOCK_SIZE); Luft 1,1 kg/m3 * 1006 J/(kg K) "
                                 "* Luftvolumen; Wand 2 mm Stahl (WAND_RHO, WAND_CP) inkl. Boden. Physikalisch, NICHT angepasst"),
        "gL": wert(p["gL"], "W/K^(4/3)", q_anp),
        "gW": wert(p["gW"], "m^2 (eps_eff * A * F)", q_anp + "; Netzteil ~0 = Anpassung findet keinen Strahlungsanteil"),
        "gLW": wert(p["gLW"], "W/K^(4/3)", q_anp),
        "sigma": wert(p["sigma"], "W/(m^2 K^4)", "Stefan-Boltzmann-Konstante (CODATA 2018)"),
        "eps_aussen": wert(p["eps_aussen"], "-", "config.EMISSIVITAET (Annahme lackiertes Blech), wie Fluent-Randbedingung"),
        "A_ab": wert(p["A_ab"], "m^2", "Schrankoberflaeche ohne Boden (Boden adiabat), aus config.B/H/T"),
        "r_int": wert(p["r_int"], "K/W", "(Tmax - Tmittel)/P aus basis_grob (Fluent)"),
        "R_matrix": wert(np.where(np.isnan(Rm), None, Rm).tolist(), "K/W",
                         "NUR INFORMATIV: Sekanten-Waermewiderstand R_ij = (T_i - T_j)/q_ij am Basis-Betriebspunkt "
                         "(35 degC, Basisleistungen); null = keine Verbindung. Simuliert wird das nichtlineare Netzwerk."),
        "G_matrix": wert(G.tolist(), "W/K", "1/R_matrix, 0 = keine Verbindung (Sekante am Basis-Betriebspunkt)"),
        "G_wand_umgebung": wert(G_umg, "W/K", "Sekante hA + Strahlung aussen am Basis-Betriebspunkt"),
        "T_ss_basis": wert(p["T_ss_basis"], "K", "stationaere ROM-Loesung Basisfall (rom.py)"),
        "A_lin": wert(p["A_lin"], "1/s", "linearisiert um T_ss_basis (rom.py), nur informativ"),
        "B_lin": wert(p["B_lin"], "K/J bzw. 1/s (Spalte T_umg)", "linearisiert, Eingaenge " + ", ".join(p["eingaenge"])),
        "stationaer": wert(stat, "P in W, T in K, hA in W/K, t in s",
                           f"Simulink: Simulation {T_STAT_H:.0f} h ab T_umg, Endwert = stationaere Loesung"),
        "transient": wert(trans, "P in W, T in K, hA in W/K, t in s", "siehe quelle im Fall"),
    }
    pj = EXPORT / "rom_parameter.json"
    pj.write_text(json.dumps(par, indent=2, ensure_ascii=True), encoding="utf-8")

    # ---------------- Pruefungen am neu geladenen JSON
    q = json.loads(pj.read_text(encoding="utf-8"))
    pruef = []

    def check(name, ok, info=""):
        pruef.append((name, ok, info))

    n2 = R.Netz(np.r_[np.log(np.maximum(np.r_[q["gL"]["wert"], q["gW"]["wert"], q["gLW"]["wert"]], 1e-300)), np.zeros(4)])
    alle = np.r_[q["gL"]["wert"], q["gW"]["wert"], q["gLW"]["wert"], q["C"]["wert"]]
    check("alle Koeffizienten endlich und >= 0", bool(np.all(np.isfinite(alle)) and np.all(alle >= 0)))
    check("6 Knoten, 4 Bauteile", len(q["knoten"]["wert"]) == 6 and len(q["gL"]["wert"]) == 4)
    dC = np.max(np.abs(np.array(q["C"]["wert"]) - R.kapazitaeten()) / R.kapazitaeten())
    check("C = rho c V (config.py)", dC < 1e-12, f"rel. Abw. {dC:.1e}")
    Tss = n2.stationaer(np.array(q["stationaer"]["wert"][0]["P_W"]), q["stationaer"]["wert"][0]["T_umg_K"], basis["h"])
    dT = np.max(np.abs(Tss - np.array(q["T_ss_basis"]["wert"])))
    check("Betriebspunkt aus JSON reproduziert", dT < 1e-6, f"max {dT:.1e} K")
    for z in q["stationaer"]["wert"]:
        check(f"hA {z['name']} plausibel (5..20 W/K)", 5 < z["hA_aussen_W_per_K"] < 20, f"{z['hA_aussen_W_per_K']:.3f} W/K")
    Gq = np.array(q["G_matrix"]["wert"])
    check("G-Matrix symmetrisch, Diagonale 0", np.allclose(Gq, Gq.T) and np.all(np.diag(Gq) == 0))

    # ---------------- Python-ROM fuer dieselben Faelle (mit hA statt h je Seite: identisch, da nur Summe eingeht)
    def rhs(T, P, T_U, hAw):
        h1 = {s: hAw / R.A_AB for s in R.SEITE_A}      # gleichmaessig verteilt -> gleiche Summe hA
        return n2.fluesse(T, P, T_U, h1) / np.array(q["C"]["wert"])

    zeilen = []
    log_stat = []
    for z in q["stationaer"]["wert"]:
        P, T_U, hAw = np.array(z["P_W"]), z["T_umg_K"], z["hA_aussen_W_per_K"]
        t1 = time.perf_counter()
        h1 = {s: hAw / R.A_AB for s in R.SEITE_A}
        T = n2.stationaer(P, T_U, h1)
        dt = time.perf_counter() - t1
        qWU = hAw * (T[5] - T_U) + q["eps_aussen"]["wert"] * R.SIGMA * R.A_AB * (T[5] ** 4 - T_U ** 4)
        check(f"Energiebilanz ROM {z['name']}", abs(qWU - P.sum()) / P.sum() < 1e-6,
              f"Q_aussen {qWU:.4f} W / P {P.sum():.1f} W")
        err = T - np.array(z["T_cfd_K"])
        log_stat.append(f"  {z['name']:14s} ({z['rolle']:11s}) max |ROM-CFD| {np.max(np.abs(err)):.2f} K  "
                        + " ".join(f"{k} {e:+.2f}" for k, e in zip(KN, err)) + f" | {1e3 * dt:.1f} ms")
        zeilen.append([z["name"], "stationaer", z["t_ende_s"]] + [round(v - 273.15, 5) for v in T] + [dt])
        # Kontrolle: Zeitintegration bis t_ende landet auf derselben stationaeren Loesung (wie in Simulink)
        iv = solve_ivp(lambda t, T: rhs(T, P, T_U, hAw), (0, z["t_ende_s"]), np.array(z["T0_K"]),
                       method="Radau", rtol=1e-8, atol=1e-8)
        check(f"Zeitintegration {T_STAT_H:.0f} h = stationaer ({z['name']})",
              np.max(np.abs(iv.y[:, -1] - T)) < 1e-4, f"{np.max(np.abs(iv.y[:, -1] - T)):.1e} K")

    tz = q["transient"]["wert"]
    tt = np.arange(0.0, tz["t_ende_s"] + 1e-9, tz["dt_ausgabe_s"])
    t1 = time.perf_counter()
    iv = solve_ivp(lambda t, T: rhs(T, np.array(tz["P_W"]), tz["T_umg_K"], tz["hA_aussen_W_per_K"]),
                   (0, tz["t_ende_s"]), np.array(tz["T0_K"]), method="Radau", t_eval=tt, rtol=1e-8, atol=1e-8)
    dt_tr = time.perf_counter() - t1
    for j, t in enumerate(iv.t):
        zeilen.append([tz["name"], "transient", t] + [round(v - 273.15, 5) for v in iv.y[:, j]] + [dt_tr])
    spalten = {"Umrichter": "tavg-umrichter", "Netzteil": "tavg-netzteil", "Schuetz": "tavg-schuetz",
               "Trafo": "tavg-trafo", "Luft": "tavg-luft"}
    dd = d.set_index("t").loc[tt]
    fehl = {k: float(np.max(np.abs(iv.y[i] - dd[c].to_numpy()))) for i, (k, c) in enumerate(spalten.items())}

    with open(EXPORT / "referenz_python.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["fall", "art", "t_s"] + [f"T_{k}_C" for k in KN] + ["rechenzeit_s"])
        wr.writerows(zeilen)

    log.append(f"rom_export erzeugt ({time.strftime('%Y-%m-%d %H:%M')}): rom_parameter.json, referenz_python.csv "
               f"({len(q['stationaer']['wert'])} stationaere Faelle + transient {len(tt)} Zeitpunkte)")
    log.append("Pruefungen JSON/Parameter:")
    for n, ok, info in pruef:
        log.append(f"  [{'OK ' if ok else 'FEHLER'}] {n} {info}")
    log.append("Python-ROM (aus JSON) gegen CFD, stationaer (Knotenmittel):")
    log += log_stat
    log.append("Python-ROM gegen CFD, transient 8 h: max |ROM-CFD| "
               + ", ".join(f"{k} {v:.2f} K" for k, v in fehl.items())
               + f" | ROM {1e3 * dt_tr:.0f} ms, CFD {tr_json['laufzeit_s'] / 60:.1f} min")
    m = pruefe_matlab_skript(EXPORT / "rom_simulink.m")
    check("rom_simulink.m: nur Bloecke der Simulink-Basisbibliothek", not m["fremd"], ", ".join(m["bloecke"]))
    check("rom_simulink.m: keine Toolbox-Funktionen", not m["toolbox"], ", ".join(m["toolbox"]))
    check("rom_simulink.m: reiner ASCII-Text (Zeichensatz egal)", m["nicht_ascii"] == 0)
    check("rom_simulink.m: Simulink-Pruefung mit exist('new_system') ohne 'file'", m["exist_ok"])
    log.append("Pruefungen rom_simulink.m (statisch):")
    for n, ok, info in pruef[-4:]:
        log.append(f"  [{'OK ' if ok else 'FEHLER'}] {n} {info}")
    log.append("  verwendete Funktionen (Einfuehrung): " + ", ".join(m["funktionen"]) + " -> alle <= R2022a")
    log.append("  matlab.codetools.requiredFilesAndProducts (R2025a, hier ohne Simulink): nur 'MATLAB'; "
               "Simulink-Aufrufe kann es ohne installiertes Simulink nicht zuordnen.")
    log.append("  MATLAB-Selbsttest rom_simulink('ohne_simulink') (gleiche Fcn-Ausdruecke, ode15s): siehe "
               "rom_export/test_ohne_simulink/rom_simulink_info.txt. Ein echter Lauf in R2022a/Simulink ist NICHT getestet.")
    log.append("Simulink-Lauf: nur auf einem Rechner mit Simulink; Auswertung und Status mit rom_vergleich.py "
               "(logs/rom_vergleich.txt). Dieser Generator fuehrt Simulink NICHT aus.")
    (HIER / "logs" / "rom_export.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))
    if not all(ok for _, ok, _ in pruef):
        raise SystemExit("Pruefung fehlgeschlagen")


if __name__ == "__main__":
    main()
