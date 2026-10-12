"""plotstil.py - gemeinsamer Stil fuer alle Diagramme (matplotlib).

Farben: feste Zuordnung je Bauteil (Identitaet folgt dem Objekt, nie dem Rang),
Temperaturfelder mit einer einfarbigen sequentiellen Skala (hell -> dunkel).
"""
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

FLAECHE = "#fcfcfb"
TEXT = "#0b0b0b"
TEXT2 = "#52514e"
GITTER = "#e4e3df"

# Kategorische Reihenfolge (validierte Referenzpalette, Slots 1-8)
KATEGORIE = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
BAUTEIL_FARBE = {"Umrichter": KATEGORIE[0], "Netzteil": KATEGORIE[1],
                 "Schuetz": KATEGORIE[2], "Trafo": KATEGORIE[3]}
NEUTRAL = "#8a8984"

# Sequentiell (Temperatur): ein Farbton, hell -> dunkel
T_CMAP = LinearSegmentedColormap.from_list(
    "temperatur", ["#fdf1e6", "#f8c9a3", "#ef9a63", "#d9622e", "#a83d16", "#6b220a"])
# Sequentiell (Geschwindigkeit): blau
V_CMAP = LinearSegmentedColormap.from_list(
    "geschw", ["#f2f6fc", "#bcd3f1", "#7aa9e6", "#2a78d6", "#1a4f95", "#0d2b55"])


def stil():
    plt.rcParams.update({
        "figure.facecolor": FLAECHE, "axes.facecolor": FLAECHE, "savefig.facecolor": FLAECHE,
        "axes.edgecolor": GITTER, "axes.labelcolor": TEXT2, "xtick.color": TEXT2, "ytick.color": TEXT2,
        "text.color": TEXT, "axes.titlecolor": TEXT, "axes.grid": True, "grid.color": GITTER,
        "grid.linewidth": 0.6, "axes.spines.top": False, "axes.spines.right": False,
        "lines.linewidth": 2.0, "font.size": 10, "axes.titlesize": 11, "legend.frameon": False,
        "figure.dpi": 150, "savefig.dpi": 150, "savefig.bbox": "tight",
    })


def endlabels(ax, eintraege, abstand, dx):
    """Mehrere Endlabels (x, y, text, farbe) ohne Ueberlappung: Text um dx rechts vom Endpunkt,
    vertikal mit Mindestabstand 'abstand' (Datenkoordinaten)."""
    eintraege = sorted(eintraege, key=lambda e: e[1])
    lagen = []
    for x, y, text, farbe in eintraege:
        yt = max(y, lagen[-1] + abstand) if lagen else y
        lagen.append(yt)
        ax.plot([x], [y], "o", ms=5, color=farbe, zorder=5)
        ax.annotate(text, (x, y), xytext=(x + dx, yt), textcoords="data", va="center", color=TEXT2,
                    fontsize=9)
    return lagen


def endlabel(ax, x, y, text, farbe):
    """Direktes Label am Linienende (Text in Textfarbe, Farbpunkt traegt die Identitaet)."""
    ax.plot([x], [y], "o", ms=5, color=farbe, zorder=5)
    ax.annotate(text, (x, y), xytext=(6, 0), textcoords="offset points", va="center",
                color=TEXT2, fontsize=9)
