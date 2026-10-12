"""
vv_laeufe.py - Phase 3: Sensitivitaets- und Modellvergleichslaeufe (Netz mittel) + Auswertung.

Aufruf:
  python vv_laeufe.py rechnen     # alle fehlenden Laeufe nacheinander (bereits fertige werden uebersprungen)
  python vv_laeufe.py auswerten   # Tabelle + Bilder aus laeufe/*/ergebnis.json
Ergebnis: ergebnisse/vv/vv_tabelle.csv, vv_sensitivitaet.png, logs/vv_zusammenfassung.txt
"""
import csv
import json
import subprocess
import sys
from pathlib import Path

import config as C

HIER = Path(__file__).resolve().parent
ZIEL = HIER / "ergebnisse" / "vv"
NETZ = "mittel"
BASIS = f"basis_{NETZ}"

# Fall -> (Beschreibung, --set Liste)
FAELLE = {
    "vv_eps06":      ("alle Emissivitaeten 0,6 (statt 0,9)", ["eps_innen=0.6", "eps_bauteil=0.6", "eps_aussen=0.6"]),
    "vv_eps_aussen03": ("nur aussen eps 0,3 (blankes/verzinktes Blech)", ["eps_aussen=0.3"]),
    "vv_k05":        ("Bauteil-k x 0,5", ["kfaktor=0.5"]),
    "vv_k2":         ("Bauteil-k x 2", ["kfaktor=2"]),
    "vv_tumg25":     ("Umgebung 25 °C", ["T_umg=25"]),
    "vv_tumg45":     ("Umgebung 45 °C", ["T_umg=45"]),
    "vv_h05":        ("Aussen-h x 0,5", ["hfaktor=0.5"]),
    "vv_h15":        ("Aussen-h x 1,5", ["hfaktor=1.5"]),
    "vv_laminar":    ("laminar statt k-omega SST", ["turb=laminar"]),
    "vv_boussinesq": ("Boussinesq statt incompressible-ideal-gas", ["dichte=boussinesq"]),
}


def fertig(fall):
    f = HIER / "laeufe" / fall / "ergebnis.json"
    return f.exists() and json.loads(f.read_text(encoding="utf-8")).get("ok")


def rechnen(auswahl=None):
    for fall, (text, setz) in FAELLE.items():
        if auswahl and fall not in auswahl:
            continue
        if fertig(fall):
            print(f"{fall}: schon fertig, uebersprungen", flush=True)
            continue
        print(f"{fall}: {text} ...", flush=True)
        with open(HIER / "logs" / f"solve_{fall}.out", "w", encoding="utf-8", errors="replace") as fh:
            r = subprocess.run([sys.executable, "-u", str(HIER / "solve_schaltschrank.py"), "--netz", NETZ,
                                "--fall", fall, "--set", *setz], cwd=HIER, stdout=fh, stderr=subprocess.STDOUT)
        print(f"{fall}: {'ok' if r.returncode == 0 else 'FEHLER (Exit ' + str(r.returncode) + ')'}", flush=True)


def auswerten():
    import plotstil as PS
    PS.stil()
    plt = PS.plt
    ZIEL.mkdir(parents=True, exist_ok=True)
    b = json.loads((HIER / "laeufe" / BASIS / "ergebnis.json").read_text(encoding="utf-8"))
    bauteile = list(C.KOMPONENTEN)
    zeilen = []
    for fall, (text, _) in ({BASIS: ("Basis (SST, ideal-gas, eps 0,9, 35 °C)", [])} | FAELLE).items():
        f = HIER / "laeufe" / fall / "ergebnis.json"
        if not f.exists():
            continue
        e = json.loads(f.read_text(encoding="utf-8"))
        if not e.get("ok"):
            zeilen.append(dict(fall=fall, beschreibung=text, status="FEHLER"))
            continue
        z = dict(fall=fall, beschreibung=text, status="ok" if e["konvergiert"] else "nicht stationaer",
                 bilanz_proz=round(e["bilanzfehler_proz"], 3), iterationen=e["iterationen"],
                 T_luft=round(e["Tluft_mittel_C"], 2), dT_luft=round(e["Tluft_mittel_C"] - b["Tluft_mittel_C"], 2),
                 T_wand=round(e["Tw_cfd_C"], 2), T_wand_hand=round(e["Tw_hand_C"], 2))
        for n in bauteile:
            z[f"Tmax_{n}"] = round(e["Tmax_C"][n], 2)
            z[f"dTmax_{n}"] = round(e["Tmax_C"][n] - b["Tmax_C"][n], 2)
        zeilen.append(z)
    felder = []
    for z in zeilen:
        for k in z:
            if k not in felder:
                felder.append(k)
    with open(ZIEL / "vv_tabelle.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=felder)
        wr.writeheader()
        wr.writerows(zeilen)
    # Bild: Aenderung der Bauteil-Tmax gegenueber Basis je Fall (gruppierte horizontale Balken)
    ok = [z for z in zeilen if z.get("status") in ("ok", "nicht stationaer") and z["fall"] != BASIS]
    if ok:
        import numpy as np
        fig, ax = plt.subplots(figsize=(7.0, 0.55 * len(ok) + 1.6))
        y = np.arange(len(ok))
        hb = 0.8 / len(bauteile)
        for i, n in enumerate(bauteile):
            ax.barh(y + (i - 1.5) * hb, [z[f"dTmax_{n}"] for z in ok], height=hb * 0.9,
                    color=PS.BAUTEIL_FARBE[n], label=n, edgecolor=PS.FLAECHE, linewidth=0.5)
        ax.set_yticks(y, [z["beschreibung"] for z in ok], fontsize=8)
        ax.axvline(0, color=PS.TEXT2, lw=0.8)
        ax.invert_yaxis()
        ax.set_xlabel("Aenderung Tmax gegenueber Basis in K")
        ax.set_title("Sensitivitaet und Modellvergleich (Netz mittel)", fontsize=9)
        ax.legend(fontsize=8, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.06 - 0.3 / len(ok)))
        ax.grid(axis="y", visible=False)
        fig.savefig(ZIEL / "vv_sensitivitaet.png")
        plt.close(fig)
    log = ["V&V-Laeufe (Netz mittel), Aenderung Tmax [K] gegenueber Basis:"]
    for z in zeilen:
        if "Tmax_Umrichter" not in z:
            log.append(f"  {z['fall']:16s} {z['status']}")
            continue
        log.append(f"  {z['fall']:16s} {z['status']:16s} Bilanz {z['bilanz_proz']:+6.2f} %  Luft {z['T_luft']:6.2f} C"
                   + "".join(f"  {n[:4]} {z[f'Tmax_{n}']:6.1f} ({z[f'dTmax_{n}']:+5.1f})" for n in bauteile))
    (HIER / "logs" / "vv_zusammenfassung.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "rechnen":
        rechnen(sys.argv[2:])
    else:
        auswerten()
