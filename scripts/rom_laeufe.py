"""
rom_laeufe.py - Phase 5: Fluent-Laeufe fuer das ROM (Netz grob, gleiche Physik wie Basis).

Stationaer (Identifikation der Leitwerte):
  basis_grob     alle Bauteile (aus Phase 2/3)
  rom_umrichter, rom_netzteil, rom_schuetz, rom_trafo : nur ein Bauteil beheizt
Stationaer (Validierung, NICHT zur Identifikation verwendet):
  rom_val        andere Leistungsverteilung und 40 degC Umgebung
Transient: rom_transient.py (Sprung aller Leistungen aus dem kalten Zustand)

Aufruf: python rom_laeufe.py
"""
import json
import subprocess
import sys
from pathlib import Path

HIER = Path(__file__).resolve().parent
NETZ = "grob"
ALLE = ("Umrichter", "Netzteil", "Schuetz", "Trafo")
FAELLE = {f"rom_{n.lower()}": [f"P_{m}=0" for m in ALLE if m != n] for n in ALLE}
FAELLE["rom_val"] = ["P_Umrichter=180", "P_Netzteil=20", "P_Schuetz=30", "P_Trafo=40", "T_umg=40"]


def main():
    for fall, setz in FAELLE.items():
        f = HIER / "laeufe" / fall / "ergebnis.json"
        if f.exists() and json.loads(f.read_text(encoding="utf-8")).get("ok"):
            print(f"{fall}: schon fertig", flush=True)
            continue
        print(f"{fall}: {setz}", flush=True)
        with open(HIER / "logs" / f"solve_{fall}.out", "w", encoding="utf-8", errors="replace") as fh:
            r = subprocess.run([sys.executable, "-u", str(HIER / "solve_schaltschrank.py"), "--netz", NETZ,
                                "--fall", fall, "--set", *setz], cwd=HIER, stdout=fh, stderr=subprocess.STDOUT)
        print(f"{fall}: {'ok' if r.returncode == 0 else 'FEHLER'}", flush=True)


if __name__ == "__main__":
    main()
