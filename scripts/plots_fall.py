"""
plots_fall.py - Bilder zu einem Fluent-Lauf aus den exportierten Daten (laeufe/<fall>/).

Aufruf:  python plots_fall.py <fall> [<fall> ...]
Ergebnis: laeufe/<fall>/bilder/*.png
  temperatur_z060.png  Temperatur + Geschwindigkeitspfeile, Ebene z = 0.06 m (durch alle Bauteile)
  temperatur_z300.png  dto. vor den Bauteilen (z = 0.30 m)
  temperatur_x200.png / temperatur_x650.png  Seitenschnitte durch Umrichter / Trafo
  schichtung.png       Lufttemperatur ueber der Hoehe (Ebene z = 0.30 m, Mittel ueber x)
  monitore.png         Verlauf der Monitorgroessen ueber die Iterationen
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle
from scipy.interpolate import griddata

import config as C
import plotstil as PS

PS.stil()
plt = PS.plt
HIER = Path(__file__).resolve().parent

# Ebene -> (horizontale Achse, Geschwindigkeitskomponente horizontal, Bezeichnung, Normale, Wert)
EBENEN = {
    "z060": ("x", "x-velocity", "Ebene z = 0,06 m (durch alle Bauteile), Blick von vorn", "z", 0.06),
    "z300": ("x", "x-velocity", "Ebene z = 0,30 m (Luftraum vor den Bauteilen), Blick von vorn", "z", 0.30),
    "x200": ("z", "z-velocity", "Ebene x = 0,20 m (durch den Umrichter), Blick von der Seite", "x", 0.20),
    "x650": ("z", "z-velocity", "Ebene x = 0,65 m (durch den Trafo), Blick von der Seite", "x", 0.65),
}


def lies_ebene(datei):
    d = pd.read_csv(datei, skipinitialspace=True)
    d.columns = [c.strip() for c in d.columns]
    return d


def bauteil_rechtecke(ax, achse, normale, wert, layout="A", em=False):
    if em:     # EM-Geometrie: Trafo in Kern/Wicklung aufgeteilt, plus Stromschienen
        import build_geometry_em as GEO
        teile = {n.replace("Trafo_", "").replace("Schiene_", ""): (o, s) for n, (o, s) in GEO.koerper().items()}
    else:
        teile = {n: (o, C.BLOCK_SIZE[n]) for n, o in C.LAYOUTS[layout].items()}
    klein = []
    for n, ((ox, oy, oz), (sx, sy, sz)) in teile.items():
        o = {"x": ox, "y": oy, "z": oz}
        g = {"x": sx, "y": sy, "z": sz}
        if not (o[normale] <= wert <= o[normale] + g[normale]):
            continue
        ax.add_patch(Rectangle((o[achse], oy), g[achse], sy, fill=False, lw=1.0,
                               ec=PS.TEXT2, ls="--", zorder=4))
        if min(g[achse], sy) < 0.03:      # sehr kleine Teile (Schienen): nur Umriss, gemeinsames Label
            klein.append((o[achse], o[achse] + g[achse], oy + sy))
            continue
        ax.text(o[achse] + g[achse] / 2, oy + sy / 2, n, ha="center", va="center",
                fontsize=8, color=PS.TEXT, zorder=6,
                bbox=dict(boxstyle="round,pad=0.2", fc=PS.FLAECHE, ec="none", alpha=0.75))
    if klein:
        ax.text((min(k[0] for k in klein) + max(k[1] for k in klein)) / 2, max(k[2] for k in klein) + 0.02,
                "Cu-Schienen L1-L3", ha="center", va="bottom", fontsize=7.5, color=PS.TEXT, zorder=6,
                bbox=dict(boxstyle="round,pad=0.2", fc=PS.FLAECHE, ec="none", alpha=0.75))


def ebene_bild(ordner, name, ziel, titel_fall, layout="A", em=False):
    achse, vkomp, titel, normale, wert = EBENEN[name]
    d = lies_ebene(ordner / f"ebene_{name}.csv")
    h = d[f"{achse}-coordinate"].to_numpy()
    y = d["y-coordinate"].to_numpy()
    T = d["temperature"].to_numpy() - 273.15
    breite = C.B if achse == "x" else C.T
    nx, ny = (320, 480) if achse == "x" else (160, 480)
    gh, gy = np.meshgrid(np.linspace(0, breite, nx), np.linspace(0, C.H, ny))
    Tg = griddata((h, y), T, (gh, gy), method="linear")
    fig, ax = plt.subplots(figsize=(6.2 if achse == "x" else 4.2, 7.0))
    lo, hi = np.nanmin(Tg), np.nanpercentile(Tg, 99.5)
    schritt = next(s for s in (0.5, 1, 2, 2.5, 5, 10) if (hi - lo) / s <= 30)
    stufen = np.arange(np.floor(lo / schritt) * schritt, np.ceil(hi / schritt) * schritt + schritt, schritt)
    cf = ax.contourf(gh, gy, Tg, levels=stufen, cmap=PS.T_CMAP, extend="max")
    cb = fig.colorbar(cf, ax=ax, shrink=0.8, pad=0.02)
    cb.set_label("Temperatur in °C")
    # Pfeile (in-plane) auf grobem Raster, nur im Luftraum (|v| > 0)
    qh, qy = np.meshgrid(np.linspace(0.02, breite - 0.02, 16 if achse == "x" else 8),
                         np.linspace(0.02, C.H - 0.02, 24))
    u = griddata((h, y), d[vkomp].to_numpy(), (qh, qy), method="linear")
    v = griddata((h, y), d["y-velocity"].to_numpy(), (qh, qy), method="linear")
    sp = np.hypot(u, v)
    maske = sp > 1e-3
    ax.quiver(qh[maske], qy[maske], u[maske], v[maske], color=PS.TEXT, alpha=0.75,
              scale=2.0, scale_units="xy", width=0.004, headwidth=3.5, zorder=5)
    bauteil_rechtecke(ax, achse, normale, wert, layout, em)
    ax.set_xlim(0, breite)
    ax.set_ylim(0, C.H)
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xlabel(f"{achse} in m")
    ax.set_ylabel("y (Hoehe) in m")
    ax.set_title(f"{titel}\n{titel_fall}\nPfeile: Geschwindigkeit in der Ebene "
                 f"(max {np.nanmax(sp):.2f} m/s)", fontsize=8.5)
    fig.savefig(ziel / f"temperatur_{name}.png")
    plt.close(fig)


def schichtung(ordner, ziel, titel_fall):
    d = lies_ebene(ordner / "ebene_z300.csv")
    d = d[(d["x-coordinate"] > 0.03) & (d["x-coordinate"] < C.B - 0.03)]   # Wandgrenzschichten aussparen
    kl = np.linspace(0, C.H, 25)
    d["klasse"] = pd.cut(d["y-coordinate"], kl)
    g = d.groupby("klasse", observed=True)["temperature"].agg(["mean", "min", "max"]) - 273.15
    ym = [iv.mid for iv in g.index]
    fig, ax = plt.subplots(figsize=(5.0, 5.0))
    ax.fill_betweenx(ym, g["min"], g["max"], color=PS.KATEGORIE[0], alpha=0.15, lw=0,
                     label="Spannweite ueber x")
    ax.plot(g["mean"], ym, color=PS.KATEGORIE[0], label="Mittel ueber x")
    ax.set_xlabel("Lufttemperatur in °C")
    ax.set_ylabel("Hoehe y in m")
    ax.set_title(f"Temperaturschichtung im Luftraum (Ebene z = 0,30 m)\n{titel_fall}", fontsize=9)
    ax.legend(loc="lower right", fontsize=8)
    fig.savefig(ziel / "schichtung.png")
    plt.close(fig)
    return float(g["mean"].iloc[-1] - g["mean"].iloc[0])


def monitore(ordner, ziel, titel_fall):
    datei = ordner / "monitor.csv"
    if not datei.exists():          # EM-Kopplung: je Kopplungsschritt eine Datei -> aneinanderhaengen
        teile = sorted(ordner.glob("monitor_k*.csv"), key=lambda f: int(f.stem.split("_k")[1]))
        if not teile:
            return
        d = pd.concat([pd.read_csv(f) for f in teile], ignore_index=True)
        d["iter"] = np.arange(1, len(d) + 1) * (d["iter"].iloc[1] - d["iter"].iloc[0])
    else:
        d = pd.read_csv(datei)
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    for n, farbe in PS.BAUTEIL_FARBE.items():
        k = f"tmax-{n.lower()}"
        if k not in d:
            continue
        ax.plot(d["iter"], d[k] - 273.15, color=farbe, label=f"Tmax {n}")
        PS.endlabel(ax, d["iter"].iloc[-1], d[k].iloc[-1] - 273.15, n, farbe)
    ax.plot(d["iter"], d["tavg-luft"] - 273.15, color=PS.NEUTRAL, label="Luft (Mittel)")
    PS.endlabel(ax, d["iter"].iloc[-1], d["tavg-luft"].iloc[-1] - 273.15, "Luft", PS.NEUTRAL)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Temperatur in °C")
    ax.set_title(f"Monitorverlauf (pseudo-transient)\n{titel_fall}", fontsize=9)
    ax.legend(fontsize=8, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.18))
    ax.set_xlim(0, d["iter"].iloc[-1] * 1.12)
    fig.savefig(ziel / "monitore.png")
    plt.close(fig)


def main():
    for fall in sys.argv[1:]:
        ordner = HIER / "laeufe" / fall
        ziel = ordner / "bilder"
        ziel.mkdir(exist_ok=True)
        e = json.loads((ordner / "ergebnis.json").read_text(encoding="utf-8"))
        p = e["parameter"]
        titel_fall = (f"{fall}: {p['turb']}, {p['dichte']}, Netz {p['netz']}, "
                      f"T_umg {p['T_umg']:.0f} °C")
        layout = p.get("layout", "A")
        if p.get("luefter"):
            titel_fall += ", Filterluefter"
        for n in EBENEN:
            if (ordner / f"ebene_{n}.csv").exists():
                ebene_bild(ordner, n, ziel, titel_fall, layout, em=str(p.get("netz", "")).startswith("em"))
        dT = schichtung(ordner, ziel, titel_fall)
        monitore(ordner, ziel, titel_fall)
        print(f"{fall}: Bilder in {ziel}  (Schichtung oben-unten {dT:+.1f} K)")


if __name__ == "__main__":
    main()
