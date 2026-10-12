"""
szenarien.py - Phase 6: Szenarienvergleich (Netz mittel).

  Basis           natuerliche Konvektion, Layout A, 35 degC   (laeufe/basis_mittel, aus Phase 2)
  45 degC         wie Basis, Umgebung 45 degC                 (laeufe/vv_tumg45, aus Phase 3)
  Layout B        Netzteil/Schuetz unten, Umrichter/Trafo oben, 35 und 45 degC
  Filterluefter   Luefter unten in der Tuer + Austrittsfilter oben, 35 und 45 degC

Aufruf:  python szenarien.py rechnen | auswerten
Ergebnis: ergebnisse/szenarien/szenarien.csv, szenarien_tmax.png, logs/szenarien.txt
"""
import csv
import json
import subprocess
import sys
from pathlib import Path

import config as C

HIER = Path(__file__).resolve().parent
ZIEL = HIER / "ergebnisse" / "szenarien"
NEU = {
    "sz_layoutB":     ["layout=B"],
    "sz_layoutB_45":  ["layout=B", "T_umg=45"],
    "sz_luefter":     ["luefter=1"],
    "sz_luefter_45":  ["luefter=1", "T_umg=45"],
}
ALLE = {
    "basis_mittel":  "Basis: natuerlich, Layout A, 35 °C",
    "vv_tumg45":     "natuerlich, Layout A, 45 °C",
    "sz_layoutB":    "natuerlich, Layout B, 35 °C",
    "sz_layoutB_45": "natuerlich, Layout B, 45 °C",
    "sz_luefter":    "Filterluefter, Layout A, 35 °C",
    "sz_luefter_45": "Filterluefter, Layout A, 45 °C",
}


def rechnen():
    for fall, setz in NEU.items():
        f = HIER / "laeufe" / fall / "ergebnis.json"
        if f.exists() and json.loads(f.read_text(encoding="utf-8")).get("ok"):
            print(f"{fall}: schon fertig", flush=True)
            continue
        print(f"{fall}: {setz}", flush=True)
        with open(HIER / "logs" / f"solve_{fall}.out", "w", encoding="utf-8", errors="replace") as fh:
            r = subprocess.run([sys.executable, "-u", str(HIER / "solve_schaltschrank.py"), "--netz", "mittel",
                                "--fall", fall, "--set", *setz], cwd=HIER, stdout=fh, stderr=subprocess.STDOUT)
        print(f"{fall}: {'ok' if r.returncode == 0 else 'FEHLER'}", flush=True)


def auswerten():
    import numpy as np
    import plotstil as PS
    PS.stil()
    plt = PS.plt
    ZIEL.mkdir(parents=True, exist_ok=True)
    bt = list(C.KOMPONENTEN)
    zeilen = []
    for fall, text in ALLE.items():
        f = HIER / "laeufe" / fall / "ergebnis.json"
        if not f.exists():
            continue
        e = json.loads(f.read_text(encoding="utf-8"))
        if not e.get("ok"):
            continue
        z = dict(fall=fall, szenario=text, stationaer=e["konvergiert"], bilanz_proz=round(e["bilanzfehler_proz"], 3),
                 T_umg=e["parameter"]["T_umg"], T_luft_mittel=round(e["Tluft_mittel_C"], 2),
                 T_luft_max=round(e["Tluft_max_C"], 2),
                 V_luefter_m3h=round(e["luefter"]["V_cfd_m3h"], 1) if "luefter" in e else "",
                 V_analytisch_m3h=round(e["luefter"]["V_analytisch_m3h"], 1) if "luefter" in e else "")
        for n in bt:
            z[f"Tmax_{n}"] = round(e["Tmax_C"][n], 2)
            z[f"dT_{n}"] = round(e["Tmax_C"][n] - e["parameter"]["T_umg"], 2)
        zeilen.append(z)
    with open(ZIEL / "szenarien.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(zeilen[0]))
        wr.writeheader()
        wr.writerows(zeilen)
    fig, ax = plt.subplots(figsize=(7.2, 0.62 * len(zeilen) + 1.8))
    y = np.arange(len(zeilen))
    hb = 0.8 / (len(bt) + 1)
    for i, n in enumerate(bt):
        ax.barh(y + (i - 2) * hb, [z[f"Tmax_{n}"] for z in zeilen], height=hb * 0.9, color=PS.BAUTEIL_FARBE[n],
                label=f"Tmax {n}", edgecolor=PS.FLAECHE, linewidth=0.5)
    ax.barh(y + 2 * hb, [z["T_luft_mittel"] for z in zeilen], height=hb * 0.9, color=PS.NEUTRAL,
            label="Luft (Mittel)", edgecolor=PS.FLAECHE, linewidth=0.5)
    for j, z in enumerate(zeilen):
        ax.text(z["Tmax_Schuetz"] + 1, j - 2 * hb + 2 * hb, f"{z['Tmax_Schuetz']:.0f}", va="center", fontsize=7,
                color=PS.TEXT2)
    ax.axvline(C.T_MAX_INNEN, color=PS.TEXT2, lw=1, ls=":")
    ax.text(C.T_MAX_INNEN + 0.5, -0.6, f"{C.T_MAX_INNEN:.0f} °C Auslegungsgrenze Innenluft (Annahme)",
            fontsize=7, color=PS.TEXT2, va="bottom")
    ax.set_yticks(y, [z["szenario"] for z in zeilen], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Temperatur in °C")
    ax.set_title("Szenarien: Bauteil-Maximaltemperaturen und mittlere Innenluft (Netz mittel)", fontsize=9)
    ax.legend(fontsize=8, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.1 - 0.25 / len(zeilen)))
    ax.grid(axis="y", visible=False)
    fig.savefig(ZIEL / "szenarien_tmax.png")
    plt.close(fig)
    log = ["Szenarien (Netz mittel):"]
    for z in zeilen:
        log.append(f"  {z['szenario']:34s} Luft {z['T_luft_mittel']:6.2f} C | " + " ".join(
            f"{n[:4]} {z[f'Tmax_{n}']:6.1f}" for n in bt) + f" | Bilanz {z['bilanz_proz']:+.2f} %"
            + (f" | V {z['V_luefter_m3h']} m3/h (analyt. {z['V_analytisch_m3h']})" if z["V_luefter_m3h"] != "" else ""))
    (HIER / "logs" / "szenarien.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))


if __name__ == "__main__":
    rechnen() if len(sys.argv) > 1 and sys.argv[1] == "rechnen" else auswerten()
