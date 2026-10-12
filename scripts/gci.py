"""
gci.py - Netzunabhaengigkeit nach dem GCI-Verfahren (Celik et al., J. Fluids Eng. 130 (2008) 078001;
Richardson-Extrapolation mit Sicherheitsfaktor 1,25 fuer drei Netze).

Groessen phi: Temperaturerhoehung ueber Umgebung (K) - Tmax je Bauteil, mittlere Luft, Wand.
Netz 1 = fein, 2 = mittel, 3 = grob; h_i = (V / N_i)^(1/3).

Aufruf:  python gci.py   (liest laeufe/basis_{grob,mittel,fein}/ergebnis.json)
Ergebnis: ergebnisse/vv/gci.csv, ergebnisse/vv/gci.png, logs/gci.txt
"""
import csv
import json
import math
from pathlib import Path

import config as C
import plotstil as PS

HIER = Path(__file__).resolve().parent
ZIEL = HIER / "ergebnisse" / "vv"
ZELLEN = {"grob": 45288, "mittel": 134913, "fein": 353186}      # aus logs/gate1_netz.txt
V = C.B * C.H * C.T


def gci3(phi1, phi2, phi3, h1, h2, h3):
    """Celik 2008, Schritte 1-5. Gibt dict mit p, phi_ext, e_a, GCI_fein, Konvergenzart."""
    r21, r32 = h2 / h1, h3 / h2
    e21, e32 = phi2 - phi1, phi3 - phi2
    out = dict(r21=r21, r32=r32, e21=e21, e32=e32)
    if abs(e21) < 1e-12 or abs(e32) < 1e-12:
        out.update(art="identisch", p=float("nan"), phi_ext=phi1, e_a=0.0, gci_fein=0.0)
        return out
    s = math.copysign(1.0, e32 / e21)
    out["art"] = "monoton" if s > 0 else "oszillierend"
    p = 2.0
    for _ in range(200):
        q = math.log((r21 ** p - s) / (r32 ** p - s))
        p_neu = abs(math.log(abs(e32 / e21)) + q) / math.log(r21)
        if abs(p_neu - p) < 1e-10:
            break
        p = p_neu
    phi_ext = (r21 ** p * phi1 - phi2) / (r21 ** p - 1)
    e_a = abs((phi1 - phi2) / phi1)
    out.update(p=p, phi_ext=phi_ext, e_a=e_a, e_ext=abs((phi_ext - phi1) / phi_ext),
               gci_fein=1.25 * e_a / (r21 ** p - 1))
    return out


def main():
    PS.stil()
    plt = PS.plt
    ZIEL.mkdir(parents=True, exist_ok=True)
    e = {n: json.loads((HIER / "laeufe" / f"basis_{n}" / "ergebnis.json").read_text(encoding="utf-8"))
         for n in ("fein", "mittel", "grob")}
    h = {n: (V / ZELLEN[n]) ** (1 / 3) for n in e}
    T_u = e["fein"]["parameter"]["T_umg"]
    groessen = {f"Tmax {b}": lambda d, b=b: d["Tmax_C"][b] - T_u for b in C.KOMPONENTEN}
    groessen["Luft Mittel"] = lambda d: d["Tluft_mittel_C"] - T_u
    groessen["Wand (Q-gew.)"] = lambda d: d["Tw_cfd_C"] - T_u
    zeilen, log = [], []
    log.append("GCI (Celik 2008), phi = Temperaturerhoehung ueber Umgebung in K")
    log.append(f"Zellen fein/mittel/grob: {ZELLEN['fein']}/{ZELLEN['mittel']}/{ZELLEN['grob']}, "
               f"h = {1e3 * h['fein']:.1f}/{1e3 * h['mittel']:.1f}/{1e3 * h['grob']:.1f} mm, "
               f"r21 = {h['mittel'] / h['fein']:.3f}, r32 = {h['grob'] / h['mittel']:.3f}")
    def ampl(d, name):
        if name.startswith("Tmax "):
            return d.get("Tmax_ampl_K", {}).get(name[5:], float("nan"))
        return d.get("Tluft_ampl_K", float("nan")) if name == "Luft Mittel" else float("nan")

    for name, fn in groessen.items():
        p1, p2, p3 = fn(e["fein"]), fn(e["mittel"]), fn(e["grob"])
        g = gci3(p1, p2, p3, h["fein"], h["mittel"], h["grob"])
        spann = max(p1, p2, p3) - min(p1, p2, p3)
        fluk = max(ampl(e[n], name) for n in e)
        # Bewertung: monoton mit plausibler Ordnung (0.5 < p < 4) -> GCI; sonst Spannweite (Celik 2008)
        gueltig = g["art"] == "monoton" and 0.5 < g["p"] < 4
        unsicher = g["gci_fein"] * p1 if gueltig else spann
        zeilen.append(dict(groesse=name, dT_fein=p1, dT_mittel=p2, dT_grob=p3, art=g["art"],
                           p=g["p"], dT_extrapoliert=g["phi_ext"], e_a_proz=100 * g["e_a"],
                           GCI_fein_proz=100 * g["gci_fein"], GCI_fein_K=g["gci_fein"] * p1,
                           Spannweite_K=spann, Schwankung_K=fluk,
                           Netzunsicherheit_K=unsicher, Methode="GCI" if gueltig else "Spannweite"))
        log.append(f"  {name:14s}: {p3:7.3f} / {p2:7.3f} / {p1:7.3f} K (grob/mittel/fein), {g['art']}, "
                   f"p = {g['p']:.2f}, GCI_fein = {100 * g['gci_fein']:.2f} % | Spannweite {spann:.3f} K, "
                   f"Schwankung +-{fluk:.3f} K -> Netzunsicherheit {unsicher:.3f} K ({'GCI' if gueltig else 'Spannweite'})")
    log.append("Bewertung: Die Unterschiede zwischen den Netzen liegen in der Groessenordnung der zeitlichen "
               "Schwankung der Monitore; die Konvergenz ist meist oszillierend und die scheinbare Ordnung p "
               "unplausibel hoch -> GCI nicht belastbar, Unsicherheit = Spannweite der drei Netze (Celik 2008).")
    with open(ZIEL / "gci.csv", "w", newline="", encoding="utf-8") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(zeilen[0]))
        wr.writeheader()
        for z in zeilen:
            wr.writerow({k: (round(v, 4) if isinstance(v, float) else v) for k, v in z.items()})
    # Bild: dT ueber h fuer die Bauteile + extrapolierter Wert bei h = 0
    fig, ax = plt.subplots(figsize=(6.4, 4.4))
    farben = dict(PS.BAUTEIL_FARBE)
    farben["Luft Mittel"] = PS.NEUTRAL
    hs = [1e3 * h[n] for n in ("fein", "mittel", "grob")]
    labels = []
    for z in zeilen:
        b = z["groesse"].replace("Tmax ", "")
        if b not in farben:
            continue
        ys = [z["dT_fein"], z["dT_mittel"], z["dT_grob"]]
        ax.plot(hs, ys, "o-", color=farben[b], ms=6)
        if z["art"] == "monoton" and math.isfinite(z["dT_extrapoliert"]):
            ax.plot([0], [z["dT_extrapoliert"]], marker="D", ms=6, color=farben[b], mfc=PS.FLAECHE)
            ax.plot([0, hs[0]], [z["dT_extrapoliert"], ys[0]], ":", color=farben[b], lw=1.2)
        labels.append((hs[-1], ys[-1], z["groesse"].replace("Tmax ", ""), farben[b]))
    PS.endlabels(ax, labels, abstand=2.8, dx=0.6)
    ax.set_xlim(-0.5, hs[-1] * 1.35)
    ax.set_xlabel("repraesentative Zellgroesse h in mm (0 = Richardson-Extrapolation)")
    ax.set_ylabel("Temperaturerhoehung ueber Umgebung in K")
    ax.set_title("Netzstudie: grob / mittel / fein (Kreise), extrapoliert (Rauten)", fontsize=9)
    fig.savefig(ZIEL / "gci.png")
    plt.close(fig)
    (HIER / "logs" / "gci.txt").write_text("\n".join(log) + "\n", encoding="utf-8")
    print("\n".join(log))


if __name__ == "__main__":
    main()
