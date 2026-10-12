"""
build_geometry_em.py - Geometrie fuer Phase 4 (EM-thermisch): schaltschrank_em.step (mm)

Gegenueber build_geometry.py:
  - Trafo-Huellquader in x aufgeteilt: Kern_L | Wicklung_L | Kern_M | Wicklung_R | Kern_R
    (Wicklungsfenster als Quader; die Vereinigung ist wieder der alte Trafo-Quader)
  - drei Cu-Stromschienen Schiene_L1..L3 (30 x 5 mm, Breitseiten zueinander)
Luftraum = Schrank minus alle Bauteile. Gemeinsame Flaechen (Kern/Wicklung) werden in Fluent
Meshing per Share Topology konform gemacht.
Aufruf: python build_geometry_em.py
"""
import cadquery as cq

import config as C

MM = 1000.0


def quader(origin, size):
    return (cq.Workplane("XY")
            .box(size[0] * MM, size[1] * MM, size[2] * MM, centered=False)
            .translate((origin[0] * MM, origin[1] * MM, origin[2] * MM)))


def koerper():
    k = {}
    for n in ("Umrichter", "Netzteil", "Schuetz"):
        k[n] = (C.BLOCK_ORIGIN[n], C.BLOCK_SIZE[n])
    ox, oy, oz = C.BLOCK_ORIGIN["Trafo"]
    sx, sy, sz = C.BLOCK_SIZE["Trafo"]
    assert abs(sum(b for _, b in C.TRAFO_TEILE_X) - sx) < 1e-9, "Trafo-Teile passen nicht in den Huellquader"
    x = ox
    for name, b in C.TRAFO_TEILE_X:
        k[f"Trafo_{name}"] = ((x, oy, oz), (b, sy, sz))
        x = round(x + b, 9)
    for ph, z0 in C.SCHIENE_Z0.items():
        k[f"Schiene_{ph}"] = ((C.SCHIENE_X[0], C.SCHIENE_Y0, z0),
                              (C.SCHIENE_X[1] - C.SCHIENE_X[0], C.SCHIENE_B, C.SCHIENE_D))
    return k


def main():
    k = koerper()
    luft = quader((0, 0, 0), (C.B, C.H, C.T))
    teile = {n: quader(o, s) for n, (o, s) in k.items()}
    for t in teile.values():
        luft = luft.cut(t)
    assy = cq.Assembly(name="Schaltschrank_EM")
    assy.add(luft, name="Luftraum")
    for n, t in teile.items():
        assy.add(t, name=n)
    assy.save("schaltschrank_em.step")
    v_luft = luft.val().Volume() / MM ** 3
    print(f"schaltschrank_em.step geschrieben: {1 + len(teile)} Koerper, Luftvolumen {v_luft * 1e3:.2f} dm3")
    for n, (o, s) in k.items():
        print(f"  {n:18s} origin {o}  size {s}")


if __name__ == "__main__":
    main()
