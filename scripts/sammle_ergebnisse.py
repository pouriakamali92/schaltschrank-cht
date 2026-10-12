"""
sammle_ergebnisse.py - Phase 7: Bilder und Tabellen fuer die Website in ergebnisse/ zusammenstellen.

Erzeugt/aktualisiert die Bilder der Einzellaeufe (plots_fall.py) und kopiert die wichtigsten nach
ergebnisse/<thema>/. Tabellen der Phasen liegen bereits in ergebnisse/vv, rom, em, szenarien.
Aufruf: python sammle_ergebnisse.py
"""
import csv
import json
import shutil
import subprocess
import sys
from pathlib import Path

import config as C

HIER = Path(__file__).resolve().parent
E = HIER / "ergebnisse"
KOPIEN = {
    "basis": ("basis_mittel", ["temperatur_z060.png", "temperatur_z300.png", "temperatur_x200.png",
                               "temperatur_x650.png", "schichtung.png", "monitore.png"]),
    "szenarien": ("sz_luefter", ["temperatur_z060.png", "temperatur_z300.png"]),
    "szenarien_layoutB": ("sz_layoutB", ["temperatur_z060.png"]),
    "em": ("em_kopplung", ["temperatur_z060.png", "temperatur_x650.png"]),
}


def main():
    faelle = sorted({f for f, _ in KOPIEN.values()} | {"basis_fein"})
    vorhanden = [f for f in faelle if (HIER / "laeufe" / f / "ergebnis.json").exists()]
    subprocess.run([sys.executable, str(HIER / "plots_fall.py"), *vorhanden], cwd=HIER, check=False)
    for thema, (fall, bilder) in KOPIEN.items():
        quelle = HIER / "laeufe" / fall / "bilder"
        ziel = E / thema.split("_")[0]
        ziel.mkdir(parents=True, exist_ok=True)
        for b in bilder:
            if (quelle / b).exists():
                name = b if thema == "basis" else f"{fall}_{b}"
                shutil.copy2(quelle / b, ziel / name)
    # Kennzahlen Basisfall als CSV
    b = json.loads((HIER / "laeufe" / "basis_mittel" / "ergebnis.json").read_text(encoding="utf-8"))
    with open(E / "basis" / "basis_kennzahlen.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.writer(fh)
        wr.writerow(["groesse", "wert", "einheit"])
        wr.writerow(["Verlustleistung gesamt", b["Q_quellen_W"], "W"])
        wr.writerow(["Energiebilanzfehler", round(b["bilanzfehler_proz"], 3), "%"])
        for n in C.KOMPONENTEN:
            wr.writerow([f"Tmax {n}", round(b["Tmax_C"][n], 2), "degC"])
            wr.writerow([f"Tmittel {n}", round(b["Tmittel_C"][n], 2), "degC"])
        wr.writerow(["Luft Mittel", round(b["Tluft_mittel_C"], 2), "degC"])
        wr.writerow(["Wand (Q-gewichtet, CFD)", round(b["Tw_cfd_C"], 2), "degC"])
        wr.writerow(["Wand Handrechnung", round(b["Tw_hand_C"], 2), "degC"])
        for s in ("boden", "decke", "links", "rechts", "rueck", "tuer"):
            wr.writerow([f"Waermestrom in Wand {s}", round(b["Q_seite_W"][s], 2), "W"])
    print("ergebnisse/ aktualisiert:")
    for p in sorted(E.rglob("*")):
        if p.is_file():
            print("  ", p.relative_to(HIER))


if __name__ == "__main__":
    main()
