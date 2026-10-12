"""
solver_info.py - liest das Netz in den Fluent-Solver und listet Zonen/Einstellungsnamen.
Nur Lesen, es wird nichts berechnet. Dient dazu, solve_schaltschrank.py exakt auf
Ihr Netz und Ihre Fluent-Version abzustimmen.

Aufruf:   D:\\Pouria\\python\\python.exe solver_info.py [grob|mittel|fein]
Ausgabe:  Konsole + solver_info.txt (bitte den GANZEN Inhalt in den Chat kopieren)
"""

import sys
from pathlib import Path

import ansys.fluent.core as pyfluent

HIER = Path(__file__).resolve().parent
stufe = sys.argv[1] if len(sys.argv) > 1 else "grob"
MESH = HIER / "mesh" / f"schaltschrank_{stufe}.msh.h5"

zeilen = []


def p(*a):
    s = " ".join(str(x) for x in a)
    print(s, flush=True)
    zeilen.append(s)


def zeige(label, fn):
    try:
        p(f"{label}: {fn()}")
    except Exception as e:
        p(f"{label}: <Fehler: {type(e).__name__}: {str(e)[:150]}>")


solver = pyfluent.launch_fluent(
    mode="solver", dimension=3, precision="double",
    processor_count=4, ui_mode="no_gui", cleanup_on_exit=True,
)
try:
    s = solver.settings
    try:
        s.file.read_mesh(file_name=MESH.as_posix())
    except Exception as e:
        p("read_mesh fehlgeschlagen, versuche read_case:", str(e)[:150])
        s.file.read_case(file_name=MESH.as_posix())
    setup = s.setup

    p("=== Zonen")
    zeige("Fluid-Zellzonen", lambda: setup.cell_zone_conditions.fluid.get_object_names())
    zeige("Solid-Zellzonen", lambda: setup.cell_zone_conditions.solid.get_object_names())
    bc = setup.boundary_conditions
    zeige("bc.child_names", lambda: bc.child_names)
    for typ in ("wall", "interface", "internal", "pressure_outlet", "pressure_inlet", "symmetry", "fan", "porous_jump"):
        zeige(f"bc.{typ}", lambda typ=typ: getattr(bc, typ).get_object_names())

    p("=== Waende: angrenzende Zelle / Shadow")
    try:
        for n in bc.wall.get_object_names():
            w = bc.wall[n]
            ac = sh = "?"
            try:
                ac = w.adjacent_cell_zone()
            except Exception:
                pass
            try:
                sh = w.shadow_face_zone()
            except Exception:
                pass
            p(f"  {n:45s} adjacent={ac}  shadow={sh}")
    except Exception as e:
        p("  Fehler:", e)

    p("=== Erlaubte Werte / Unterstrukturen")
    zeige("models.child_names", lambda: setup.models.child_names)
    zeige("models.potential", lambda: setup.models.potential.child_names if hasattr(setup.models, "potential") else "nicht vorhanden")
    zeige("models.joule/mhd", lambda: [n for n in setup.models.child_names if "mhd" in n or "joule" in n or "potential" in n or "electric" in n])
    zeige("viscous.model", lambda: setup.models.viscous.model.allowed_values())
    zeige("radiation.model", lambda: setup.models.radiation.model.allowed_values())
    zeige("general.child_names", lambda: setup.general.child_names)
    zeige("general.solver.child_names", lambda: setup.general.solver.child_names)
    zeige("operating_conditions", lambda: setup.general.operating_conditions.child_names)
    zeige("materials.child_names", lambda: setup.materials.child_names)
    zeige("materials.fluid", lambda: setup.materials.fluid.get_object_names())
    zeige("materials.solid", lambda: setup.materials.solid.get_object_names())
    zeige("air.child_names", lambda: setup.materials.fluid["air"].child_names)
    zeige("air.density.option", lambda: setup.materials.fluid["air"].density.option.allowed_values())
    zeige("air.viscosity.option", lambda: setup.materials.fluid["air"].viscosity.option.allowed_values())
    zeige("solid-mat.child_names", lambda: setup.materials.solid[setup.materials.solid.get_object_names()[0]].child_names)
    zeige("s2s", lambda: setup.models.radiation.s2s.child_names)
    zeige("viewfactor_settings", lambda: setup.models.radiation.s2s.viewfactor_settings.child_names)

    try:
        namen = bc.wall.get_object_names()
        w0 = bc.wall[namen[0]]
        zeige("wall.thermal.child_names", lambda: w0.thermal.child_names)
        zeige("wall.thermal.thermal_bc", lambda: w0.thermal.thermal_bc.allowed_values())
        zeige("wall.radiation.child_names", lambda: w0.radiation.child_names)
        zeige("wall.momentum.child_names", lambda: w0.momentum.child_names)
    except Exception as e:
        p("Wand-Info Fehler:", e)

    try:
        sol = setup.cell_zone_conditions.solid.get_object_names()
        if sol:
            z = setup.cell_zone_conditions.solid[sol[0]]
            zeige("solid-zone.child_names", lambda: z.child_names)
            zeige("solid-zone.sources.child_names", lambda: z.sources.child_names)
            zeige("solid-zone.sources.terms.child_names", lambda: z.sources.terms.child_names)
    except Exception as e:
        p("Solid-Info Fehler:", e)

    zeige("solution.child_names", lambda: s.solution.child_names)
    zeige("methods.child_names", lambda: s.solution.methods.child_names)
    zeige("methods.p_v_coupling", lambda: s.solution.methods.p_v_coupling.child_names)
    zeige("methods.discretization_scheme", lambda: s.solution.methods.discretization_scheme.child_names)
    zeige("controls.child_names", lambda: s.solution.controls.child_names)
    zeige("initialization.child_names", lambda: s.solution.initialization.child_names)
    zeige("run_calculation.child_names", lambda: s.solution.run_calculation.child_names)
    zeige("monitor.child_names", lambda: s.solution.monitor.child_names)
    zeige("reports.child_names", lambda: s.solution.report_definitions.child_names)
    zeige("pyfluent", lambda: pyfluent.__version__)
finally:
    try:
        solver.exit()
    except Exception:
        pass
    (HIER / "solver_info.txt").write_text("\n".join(zeilen), encoding="utf-8")
    print("\nGeschrieben: solver_info.txt")
