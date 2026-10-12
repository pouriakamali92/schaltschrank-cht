"""
rom_vergleich.py - Auswertung des Simulink-ROM-Laufs (Laptop) gegen Python-ROM und CFD.

Liest  rom_export/ergebnis/rom_simulink_ergebnis.csv + rom_simulink_info.txt (vom Laptop zurueckkopiert),
       rom_export/referenz_python.csv, CFD: stationaer T_cfd_K aus rom_export/rom_parameter.json,
       transient laeufe/rom_transient/transient.csv; CFD-Laufzeiten aus laeufe/*/ergebnis.json.
Schreibt ergebnisse/rom/rom_simulink_vergleich.csv und .md, logs/rom_vergleich.txt und ersetzt die
       Simulink-Bloecke (Marker <!-- SIMULINK...:START/END -->) in docs/status.md und ergebnisse/README.md.
Sperre: Nur wenn rom_simulink_info.txt "Modus: simulink" meldet, wird Simulink als ausgefuehrt eingetragen.

Aufruf: python rom_vergleich.py          (nach dem Laptop-Lauf)
        python rom_vergleich.py --test   (Probe mit rom_export/test_ohne_simulink; schreibt NUR
                                          logs/rom_vergleich_test.txt, aendert keine Doku)
"""
import json
import re
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HIER = Path(__file__).resolve().parent
EXP = HIER / "rom_export"
GRENZE_K = 0.05          # Simulink und Python rechnen dieselben Gleichungen -> Abweichung nur Solver-Toleranz
KN = ["Umrichter", "Netzteil", "Schuetz", "Trafo", "Luft", "Wand"]
TR_SPALTE = {"Umrichter": "tavg-umrichter", "Netzteil": "tavg-netzteil", "Schuetz": "tavg-schuetz",
             "Trafo": "tavg-trafo", "Luft": "tavg-luft"}


def ersetze_block(datei, marke, neu):
    txt = datei.read_text(encoding="utf-8")
    muster = re.compile(rf"(<!-- {marke}:START -->\n).*?(<!-- {marke}:END -->)", re.S)
    if not muster.search(txt):
        raise SystemExit(f"Marker {marke} fehlt in {datei}")
    datei.write_text(muster.sub(lambda m: m.group(1) + neu.rstrip() + "\n" + m.group(2), txt), encoding="utf-8")


def main():
    test = "--test" in sys.argv
    quelle = EXP / ("test_ohne_simulink" if test else "ergebnis")
    erg_csv = quelle / "rom_simulink_ergebnis.csv"
    if not erg_csv.exists():
        raise SystemExit(f"{erg_csv} fehlt - Simulink-Lauf auf dem Laptop noch nicht zurueckkopiert. "
                         "Stand bleibt: Simulink-Modell vorbereitet, noch nicht ausgefuehrt.")
    info = (quelle / "rom_simulink_info.txt").read_text(encoding="utf-8").splitlines()
    ist_simulink = any(z.startswith("Modus: simulink") for z in info)
    if not test and not ist_simulink:
        raise SystemExit("rom_simulink_info.txt meldet keinen Simulink-Lauf (Modus) - nichts eingetragen.")
    sim = pd.read_csv(erg_csv)
    py = pd.read_csv(EXP / "referenz_python.csv")
    par = json.loads((EXP / "rom_parameter.json").read_text(encoding="utf-8"))
    if len(sim) != len(py) or list(sim.fall) != list(py.fall):
        raise SystemExit("Zeilen von Simulink-Ergebnis und Python-Referenz passen nicht zusammen.")
    sp = [f"T_{k}_C" for k in KN]

    zeilen = []
    for z in par["stationaer"]["wert"] + [par["transient"]["wert"]]:
        f = z["name"]
        s, p = sim[sim.fall == f], py[py.fall == f]
        e = json.loads((HIER / "laeufe" / f / "ergebnis.json").read_text(encoding="utf-8"))
        if f == "rom_transient":
            d = pd.read_csv(HIER / "laeufe" / f / "transient.csv").set_index("t").loc[s.t_s.to_numpy()]
            cfd = np.column_stack([d[TR_SPALTE[k]].to_numpy() - 273.15 for k in KN[:5]])
            n = 5                       # Wand-Mittel im transienten Lauf nicht als eine Groesse protokolliert
        else:
            cfd = np.array(z["T_cfd_K"])[None, :] - 273.15
            n = 6
        S, P = s[sp].to_numpy(), p[sp].to_numpy()
        zeilen.append(dict(
            fall=f, art=s.art.iloc[0], rolle=z["rolle"],
            max_sim_py_K=float(np.max(np.abs(S - P))),
            max_sim_cfd_K=float(np.max(np.abs(S[:, :n] - cfd))),
            max_py_cfd_K=float(np.max(np.abs(P[:, :n] - cfd))),
            rechenzeit_simulink_s=float(s.rechenzeit_s.iloc[0]),
            rechenzeit_python_s=float(p.rechenzeit_s.iloc[0]),
            rechenzeit_cfd_min=e["laufzeit_s"] / 60))
    tab = pd.DataFrame(zeilen)
    dmax = tab.max_sim_py_K.max()
    ok = dmax < GRENZE_K
    name_sim = "Simulink" if ist_simulink else "MATLAB-Selbsttest (ohne Simulink)"
    version = next((z for z in info if z.startswith("Simulink ")), next((z for z in info if z.startswith("MATLAB")), ""))

    md = [f"| Fall | Rolle | max. Abw. {name_sim} - Python [K] | max. Abw. {name_sim} - CFD [K] | "
          "max. Abw. Python - CFD [K] | Rechenzeit Simulink [s] | Python [s] | CFD [min] |",
          "|---|---|---|---|---|---|---|---|"]
    for z in zeilen:
        md.append(f"| {z['fall']} ({z['art']}) | {z['rolle']} | {z['max_sim_py_K']:.4f} | {z['max_sim_cfd_K']:.2f} | "
                  f"{z['max_py_cfd_K']:.2f} | {z['rechenzeit_simulink_s']:.3f} | {z['rechenzeit_python_s']:.4f} | "
                  f"{z['rechenzeit_cfd_min']:.1f} |")
    kopf = (f"{name_sim}-ROM ({version}) gegen Python-ROM: max. {dmax:.4f} K "
            f"({'bestanden' if ok else 'NICHT bestanden'}, Grenze {GRENZE_K} K). Simulink-Rechenzeit inkl. "
            f"Modell-Kompilierung. Ausgewertet {time.strftime('%Y-%m-%d %H:%M')} mit rom_vergleich.py.")
    log = [kopf] + md

    if test:
        (HIER / "logs" / "rom_vergleich_test.txt").write_text(
            "PROBELAUF mit test_ohne_simulink - KEIN Simulink-Ergebnis\n" + "\n".join(log) + "\n", encoding="utf-8")
        print("\n".join(log))
        return

    ziel = HIER / "ergebnisse" / "rom"
    tab.to_csv(ziel / "rom_simulink_vergleich.csv", index=False, float_format="%.5g")
    (ziel / "rom_simulink_vergleich.md").write_text("\n".join(log) + "\n", encoding="utf-8")
    (HIER / "logs" / "rom_vergleich.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    for png in ["rom_simulink.png"]:
        (ziel / png).write_bytes((quelle / png).read_bytes())

    tr = tab[tab.art == "transient"].iloc[0]
    st_val = tab[tab.rolle == "Validierung"].iloc[0]
    st = tab[tab.art == "stationaer"]
    fit_max = st[st.rolle == "Anpassung"].max_sim_cfd_K.max()
    t_sim = (st.rechenzeit_simulink_s.min(), st.rechenzeit_simulink_s.max())
    t_cfd = (st.rechenzeit_cfd_min.min(), st.rechenzeit_cfd_min.max())
    faktor = (60 * st.rechenzeit_cfd_min / st.rechenzeit_simulink_s)
    f_tr = 60 * tr.rechenzeit_cfd_min / tr.rechenzeit_simulink_s
    laufzeit = (f"Rechenzeit Simulink {t_sim[0]:.2f}-{t_sim[1]:.2f} s je stationaerem Fall, {tr.rechenzeit_simulink_s:.2f} s "
                f"transient 8 h (inkl. Modell-Kompilierung) gegen Fluent {t_cfd[0]:.1f}-{t_cfd[1]:.1f} min stationaer und "
                f"{tr.rechenzeit_cfd_min:.1f} min transient (12 Kerne, Netz grob; laufzeit_s aus laeufe/*/ergebnis.json, "
                f"auch logs/rom.txt) -> Faktor {faktor.min():.0f}-{faktor.max():.0f} bzw. {f_tr:.0f}.")
    if ok:
        status = (f"- Simulink-ROM lief erfolgreich (Laptop, MATLAB R2024a / {version}). Modell aus Fcn- und "
                  "Integrator-Bloecken mit ode15s (KEINE Simscape-Thermal-Bloecke), per Skript erzeugt.\n"
                  f"- max |Simulink - Python-ROM| = {dmax:.4f} K: zeigt NUR die Konsistenz beider Implementierungen bei "
                  "gleichen Gleichungen und Parametern, KEINE Validierung des ROM.\n"
                  f"- ROM-Fehler gegen CFD unveraendert: Validierung {st_val.max_sim_cfd_K:.2f} K, Anpassung bis "
                  f"{fit_max:.2f} K, transient {tr.max_sim_cfd_K:.2f} K.\n"
                  f"- {laufzeit}\n- Tabelle: ergebnisse/rom/rom_simulink_vergleich.md")
        de = (f"Simulink-Modell (per Skript aus Fcn-/Integrator-Blöcken erzeugt, ode15s, keine Simscape-Blöcke) lief "
              f"erfolgreich und stimmt mit dem Python-ROM bis {dmax:.4f} K überein – das belegt nur die Konsistenz der "
              f"beiden Implementierungen, keine Validierung. Der ROM-Fehler gegen CFD bleibt {st_val.max_sim_cfd_K:.2f} K "
              f"(Validierung) bzw. bis {fit_max:.1f} K (Anpassung). Rechenzeit ca. {t_sim[1]:.1f} s je Fall gegen "
              f"{t_cfd[0]:.0f}–{t_cfd[1]:.0f} min Fluent. → `rom/rom_simulink_vergleich.md`")
        en = (f"The Simulink model (built by script from Fcn/Integrator blocks, ode15s, no Simscape blocks) ran "
              f"successfully and matches the Python ROM within {dmax:.4f} K – this only shows that both implementations "
              f"are consistent, it is not a validation. The ROM error vs. CFD remains {st_val.max_sim_cfd_K:.2f} K "
              f"(validation) and up to {fit_max:.1f} K (fit cases). Run time about {t_sim[1]:.1f} s per case vs. "
              f"{t_cfd[0]:.0f}–{t_cfd[1]:.0f} min in Fluent. → `rom/rom_simulink_vergleich.md`")
        log.append("")
        log.append("Einordnung: Simulink = Python zeigt nur Konsistenz (gleiche Gleichungen/Parameter), keine Validierung. "
                   f"ROM gegen CFD: Validierung {st_val.max_sim_cfd_K:.2f} K, Anpassung bis {fit_max:.2f} K.")
        log.append(laufzeit)
        (ziel / "rom_simulink_vergleich.md").write_text("\n".join(log) + "\n", encoding="utf-8")
        (HIER / "logs" / "rom_vergleich.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    else:
        status = (f"- Simulink-ROM ausgefuehrt, aber Abweichung zum Python-ROM {dmax:.3f} K > {GRENZE_K} K -> "
                  "Ursache klaeren (logs/rom_vergleich.txt). Nicht als erledigt werten.")
        de = "Simulink-Modell ausgeführt, Abgleich mit dem Python-ROM noch offen (Abweichung wird geklärt)."
        en = "Simulink model executed; agreement with the Python ROM still under investigation."
    ersetze_block(HIER / "docs" / "status.md", "SIMULINK", status + "\n\n" + "\n".join(md))
    ersetze_block(HIER / "ergebnisse" / "README.md", "SIMULINK-DE", de)
    ersetze_block(HIER / "ergebnisse" / "README.md", "SIMULINK-EN", en)
    print("\n".join(log))


if __name__ == "__main__":
    main()
