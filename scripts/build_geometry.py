"""
build_geometry.py - erzeugt die Schaltschrank-Geometrie als STEP-Datei (Einheit mm).

Koerper:
  - Luftraum   : Quader B x H x T MINUS Bauteile (und ggf. Oeffnungskaesten)
  - Umrichter, Netzteil, Schuetz, Trafo : massive Quader, passen exakt in die Hohlraeume
  - optional (--luefter): Oeffnung_Ein / Oeffnung_Aus = 10 mm tiefe Luftkaesten an der Tuer
    (Filterluefter unten, Austrittsfilter oben). Ihre Aussenflaeche wird eine eigene Randzone;
    im Solver werden die Kaesten mit dem Luftraum zu EINER Fluid-Zone verschmolzen.

Alle Masse aus config.py (Meter).
Aufruf: python build_geometry.py [A|B] [--luefter]
   A ohne --luefter -> schaltschrank.step (Basisgeometrie)
   sonst            -> schaltschrank_layoutB.step, schaltschrank_luefter.step, ...
"""
import sys

import cadquery as cq

from config import AUSTRITT_MITTE, B, BLOCK_SIZE, H, LAYOUTS, LUEFTER_MITTE, LUEFTER_OEFFNUNG, T

MM = 1000.0
KASTEN_TIEFE = 0.010


def quader(origin, size):
    return (
        cq.Workplane("XY")
        .box(size[0] * MM, size[1] * MM, size[2] * MM, centered=False)
        .translate((origin[0] * MM, origin[1] * MM, origin[2] * MM))
    )


def dateiname(layout="A", luefter=False):
    return "schaltschrank" + ("" if layout == "A" else f"_layout{layout}") + ("_luefter" if luefter else "") + ".step"


def main():
    layout = next((a for a in sys.argv[1:] if not a.startswith("--")), "A")
    luefter = "--luefter" in sys.argv
    origin = LAYOUTS[layout]
    teile = {n: quader(origin[n], BLOCK_SIZE[n]) for n in origin}
    if luefter:
        a = LUEFTER_OEFFNUNG / 2
        for name, (xm, ym) in (("Oeffnung_Ein", LUEFTER_MITTE), ("Oeffnung_Aus", AUSTRITT_MITTE)):
            teile[name] = quader((xm - a, ym - a, T - KASTEN_TIEFE), (LUEFTER_OEFFNUNG, LUEFTER_OEFFNUNG, KASTEN_TIEFE))
    luftraum = quader((0.0, 0.0, 0.0), (B, H, T))
    for t in teile.values():
        luftraum = luftraum.cut(t)
    assy = cq.Assembly(name="Schaltschrank")
    assy.add(luftraum, name="Luftraum")
    for name, t in teile.items():
        assy.add(t, name=name)
    datei = dateiname(layout, luefter)
    assy.save(datei)
    print(f"{datei} geschrieben (Einheit: mm), Layout {layout}, Oeffnungen: {'ja' if luefter else 'nein'}")


if __name__ == "__main__":
    main()
