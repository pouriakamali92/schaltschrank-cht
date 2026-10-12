"""
pruefe_netz.py - prueft automatisch, ob das Netz in Fluent die richtigen Zonen hat.
Erwartung: Fluid = ['luftraum-luftraum'], Solid = 4 Bauteile.
Aufruf:  python pruefe_netz.py [grob|mittel|fein]    (Exit-Code 0 = ok, 1 = Fehler)
"""
import sys
from pathlib import Path

import ansys.fluent.core as pyfluent

HIER = Path(__file__).resolve().parent
stufe = sys.argv[1] if len(sys.argv) > 1 else "grob"
MESH = HIER / "mesh" / f"schaltschrank_{stufe}.msh.h5"
BAUTEILE = {"umrichter-umrichter", "netzteil-netzteil", "schuetz-schuetz", "trafo-trafo"}

ok = False
solver = pyfluent.launch_fluent(
    mode="solver", dimension=3, precision="double",
    processor_count=4, ui_mode="no_gui", cleanup_on_exit=True, cwd=str(HIER / "mesh"),
)
try:
    s = solver.settings
    try:
        s.file.read_mesh(file_name=MESH.as_posix())
    except Exception:
        s.file.read_case(file_name=MESH.as_posix())
    czc = s.setup.cell_zone_conditions
    fluid = sorted(czc.fluid.get_object_names())
    solid = sorted(czc.solid.get_object_names())
    print(f"ZONEN fluid={fluid} solid={solid}")
    ok = fluid == ["luftraum-luftraum"] and set(solid) == BAUTEILE
    print("ZONENPRUEFUNG:", "OK" if ok else "FEHLER (Luftraum ist nicht die einzige Fluid-Zone)")
finally:
    try:
        solver.exit()
    except Exception:
        pass
sys.exit(0 if ok else 1)
