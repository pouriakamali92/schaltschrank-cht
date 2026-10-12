"""audit_fluent_lesen.py - nur LESEN: Groessen aus vorhandenen Fluent-Laeufen fuer plausibilitaet_audit.py.
Liest laeufe/basis_mittel und laeufe/sz_luefter (cas/dat), legt Report-Definitionen nur in der Sitzung an
(Praefix aud-), speichert NICHTS in den Case/Data-Dateien. Ausgabe: logs/audit_fluent.json.
Aufruf: python audit_fluent_lesen.py   (Fluent, 4 Kerne, ca. 3 min)"""
import json
import tempfile
from pathlib import Path

import ansys.fluent.core as pyfluent

PROJ = Path(r"D:\Pouria\projekts\schaltschrank")
AUS = PROJ / "logs" / "audit_fluent.json"
erg = {}


def compute(s, defs):
    roh = s.settings.solution.report_definitions.compute(report_defs=defs)
    out = {}
    for e in roh:
        for k, v in e.items():
            out[k.replace("aud-", "", 1)] = v[0] if isinstance(v, (list, tuple)) else v
    return out


def surf(s, name, typ, feld, flaechen):
    rd = s.settings.solution.report_definitions.surface
    name = "aud-" + name
    rd.create(name=name)
    rd[name].report_type = typ
    if feld:
        rd[name].field = feld
    rd[name].surface_names = flaechen
    return name


def flux(s, name, typ, raender):
    rd = s.settings.solution.report_definitions.flux
    name = "aud-" + name
    rd.create(name=name)
    rd[name].report_type = typ
    rd[name].boundaries = raender
    return name


s = pyfluent.launch_fluent(mode="solver", dimension=3, precision="double", processor_count=4, ui_mode="no_gui",
                           cleanup_on_exit=True, cwd=tempfile.mkdtemp(prefix="audit_fluent_"))
try:
    # --- Basis mittel: Bauteiloberflaechen
    f = PROJ / "laeufe" / "basis_mittel" / "basis_mittel.cas.h5"
    s.settings.file.read_case_data(file_name=f.as_posix())
    waende = list(s.settings.setup.boundary_conditions.wall.keys())
    erg["waende"] = waende
    b = {}
    for n in ["umrichter", "netzteil", "schuetz", "trafo"]:
        z = [w for w in waende if w.startswith(f"luftraum-luftraum-{n}") and not w.endswith("shadow")]
        defs = [surf(s, f"a-{n}", "surface-area", None, z),
                surf(s, f"t-{n}", "surface-areaavg", "temperature", z),
                flux(s, f"qt-{n}", "flux-heattransfer", z),
                flux(s, f"qr-{n}", "flux-radheattransfer", z)]
        b[n] = dict(zonen=z, **compute(s, defs))
    erg["basis_mittel"] = b
    # --- Luefter: Massenstrom und Temperatur an den Oeffnungen
    f = PROJ / "laeufe" / "sz_luefter" / "sz_luefter.cas.h5"
    s.settings.file.read_case_data(file_name=f.as_posix())
    bc = s.settings.setup.boundary_conditions
    alle = []
    for typ in ("intake_fan", "outlet_vent"):
        try:
            alle += list(getattr(bc, typ).keys())
        except Exception as e:  # noqa
            erg[f"fehler_{typ}"] = str(e)
    erg["oeffnungen"] = alle
    ein = [x for x in alle if "ein" in x]
    aus = [x for x in alle if "aus" in x]
    defs = [flux(s, "m-ein", "flux-massflow", ein), flux(s, "m-aus", "flux-massflow", aus),
            flux(s, "q-ein", "flux-heattransfer", ein), flux(s, "q-aus", "flux-heattransfer", aus),
            surf(s, "tm-ein", "surface-massavg", "temperature", ein),
            surf(s, "tm-aus", "surface-massavg", "temperature", aus),
            surf(s, "ta-aus", "surface-areaavg", "temperature", aus)]
    erg["sz_luefter"] = compute(s, defs)
    erg["sz_luefter"]["cp_luft"] = s.settings.setup.materials.fluid["air"].specific_heat.value()
finally:
    AUS.write_text(json.dumps(erg, indent=2, default=str), encoding="utf-8")
    s.exit()
print(json.dumps(erg, indent=1, default=str)[:4000])
