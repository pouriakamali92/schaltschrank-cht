"""
mesh_schaltschrank.py - Phase 3: Netz fuer den Schaltschrank mit Fluent Meshing (PyFluent).

Ablauf (Watertight-Workflow, Fluid + 4 Festkoerper):
  Geometrie importieren -> Oberflaechennetz -> Geometrie beschreiben (Fluid + Solid)
  -> Share Topology -> Randbedingungen/Regionen -> Volumennetz (Poly-Hexcore) -> speichern

Aufruf:   D:\\Pouria\\python\\python.exe mesh_schaltschrank.py [grob|mittel|fein]
Eingabe:  schaltschrank.step im selben Ordner (Einheit mm)
Ausgabe:  mesh\\schaltschrank_<stufe>.msh.h5

STATUS: mit Fluent 2024 R2 / PyFluent 0.42.2 getestet (alle Stufen + EM-, Layout-B- und
Luefter-Geometrie, Gate 1 erfuellt). Fallstricke (Regionstyp, Groessen in mm, set_state,
Ueberschreiben der Netzdatei) siehe docs/log.md, Phase 1.
Optional: Umgebungsvariable MESH_OVERRIDE (JSON) fuer andere Geometrie/Groessen, siehe main().

Enthalten: Prismenschichten (schlagen sie fehl, wird ohne Prismen weitergemacht und es steht
im Log) und Qualitaetsverbesserung (Skewness/Orthogonal Quality im Wechsel) mit Endkontrolle.
"""

import json
import os
import sys
import time
import traceback
from pathlib import Path

import ansys.fluent.core as pyfluent

HIER = Path(__file__).resolve().parent
STEP = HIER / "schaltschrank.step"
MESH_DIR = HIER / "mesh"
NPROC = 4

# Groessen in MM: Die Workflow-Groessen sind in der Geometrieeinheit (mm), nicht in mesh_unit
# (Vorgabewerte von Fluent: max 30, min 1.17). Verfeinerungsfaktor je Stufe ca. 1.4 (fuer GCI).
STUFEN = {
    "grob":   dict(min_size=5.0, max_size=40.0, hex_max=40.0),
    "mittel": dict(min_size=2.8, max_size=22.0, hex_max=22.0),
    "fein":   dict(min_size=1.8, max_size=14.0, hex_max=14.0),
}
OQ_GRENZE = 0.15
WACHSTUM = 1.2
BL_SCHICHTEN = 7           # Anzahl Prismenschichten
BL_WACHSTUM = 1.2          # Wachstumsrate der Prismenschichten
BL_UEBERGANG = 0.272       # Transition Ratio (smooth-transition, Fluent-Standard)
SKEWNESS_GRENZE = 0.85     # Zellen mit hoeherer Skewness werden im Verbesserungsschritt bearbeitet
ZELLEN_PRO_SPALT = 3        # Spalt Bauteil-Rueckwand = 0,02 m
KRUEMMUNGSWINKEL = 18
SOLID_ERKENNUNG = "The geometry consists of both fluid and solid regions and/or voids"
FLUID = ["luftraum-luftraum"]   # Fluid-Regionen (Luefter-Variante: + Oeffnungskaesten)


def schritt(name, fn):
    t0 = time.time()
    print(f"\n=== {name} ...", flush=True)
    try:
        ergebnis = fn()
    except Exception:
        print(f"!!! Fehler in Schritt: {name}", flush=True)
        traceback.print_exc()
        raise
    print(f"    ok ({time.time() - t0:.1f} s)", flush=True)
    return ergebnis


def task(wt, *namen):
    """Erste vorhandene Aufgabe aus den Kandidatennamen holen."""
    for n in namen:
        try:
            return getattr(wt, n)
        except Exception:
            continue
    raise RuntimeError(
        f"Aufgabe nicht gefunden: {namen}. Vorhandene Aufgaben: {wt.task_names()}"
    )


def setze(obj, wert, *namen):
    """Setzt das erste vorhandene Attribut aus den Kandidatennamen."""
    fehler = []
    for n in namen:
        try:
            setattr(obj, n, wert)
            return n
        except Exception as e:
            fehler.append(f"{n}: {str(e)[:80]}")
    raise RuntimeError(f"Kein Attribut setzbar {namen} = {wert!r}: {fehler}")


def main():
    global ZELLEN_PRO_SPALT, BL_SCHICHTEN, STEP, FLUID
    stufe = sys.argv[1] if len(sys.argv) > 1 else "mittel"
    if stufe not in STUFEN:
        sys.exit(f"Unbekannte Stufe '{stufe}'. Erlaubt: {', '.join(STUFEN)}")
    # Optionale Uebersteuerung (Netzstudien, andere Geometrie), z.B.
    # MESH_OVERRIDE='{"min_size":3,"max_size":24,"gap":4,"layers":6,"suffix":"_t1","step":"x.step","name":"em"}'
    ov = json.loads(os.environ.get("MESH_OVERRIDE", "{}") or "{}")
    s = dict(STUFEN[stufe])
    for k in ("min_size", "max_size", "hex_max"):
        if k in ov:
            s[k] = float(ov[k])
    if "max_size" in ov and "hex_max" not in ov:
        s["hex_max"] = s["max_size"]
    ZELLEN_PRO_SPALT = ov.get("gap", ZELLEN_PRO_SPALT)
    FLUID = sorted(ov.get("fluid", FLUID))
    BL_SCHICHTEN = int(ov.get("layers", BL_SCHICHTEN))
    if "step" in ov:
        STEP = (HIER / ov["step"]).resolve()
    if not STEP.exists():
        sys.exit(f"STEP-Datei nicht gefunden: {STEP}")
    MESH_DIR.mkdir(exist_ok=True)
    mesh_datei = MESH_DIR / f"{ov.get('name', 'schaltschrank')}_{stufe}{ov.get('suffix', '')}.msh.h5"
    print(f"Netzstufe: {stufe}  {s}  Spalt-Zellen {ZELLEN_PRO_SPALT}  Prismen {BL_SCHICHTEN}  STEP {STEP.name}")

    meshing = schritt(
        "Fluent Meshing starten",
        lambda: pyfluent.launch_fluent(
            mode="meshing",
            dimension=3,
            precision="double",
            processor_count=NPROC,
            ui_mode="no_gui",
            cleanup_on_exit=True,
            cwd=str(MESH_DIR),          # Fluent-Transkripte (.trn) nicht im Projektordner ablegen
        ),
    )

    try:
        wt = schritt("Watertight-Workflow initialisieren", lambda: meshing.watertight())
        print("    Aufgaben:", wt.task_names())

        # 1 --- Import
        def importieren():
            t = task(wt, "import_geometry")
            t.file_name = str(STEP)
            t.length_unit = "mm"      # Einheit im STEP
            t.mesh_unit = "m"         # Netz/Groessen in m
            try:
                t.cad_import_options.one_object_per = "body"
            except Exception as e:
                print("    Hinweis: one_object_per nicht gesetzt:", e)
            t()
        schritt("Geometrie importieren", importieren)
        print("    Aufgaben:", wt.task_names())

        # 2 --- Oberflaechennetz
        def oberflaeche():
            t = task(wt, "create_surface_mesh", "generate_the_surface_mesh_wtm")
            # Einzelne Attribut-Zuweisungen wurden in 2024 R2 NICHT uebernommen -> set_state
            t.arguments.set_state({"cfd_surface_mesh_controls": {
                "min_size": s["min_size"], "max_size": s["max_size"], "growth_rate": WACHSTUM,
                "curvature_normal_angle": KRUEMMUNGSWINKEL, "cells_per_gap": ZELLEN_PRO_SPALT}})
            c = t.arguments().get("cfd_surface_mesh_controls", {})
            print("    Oberflaechennetz-Parameter:", {k: c.get(k) for k in
                  ("min_size", "max_size", "growth_rate", "cells_per_gap", "curvature_normal_angle")})
            if abs(c.get("max_size", 0) - s["max_size"]) > 1e-9:
                raise RuntimeError("max_size wurde nicht uebernommen")
            t()
        schritt("Oberflaechennetz erzeugen", oberflaeche)

        # 3 --- Geometrie beschreiben: Fluid + Solid
        def beschreiben():
            d = task(wt, "describe_geometry")
            d.update_child_tasks(setup_type_changed=False)
            try:
                print("    Erlaubte setup_type-Werte:", d.setup_type.allowed_values())
            except Exception:
                pass
            d.setup_type = SOLID_ERKENNUNG
            d.update_child_tasks(setup_type_changed=True)
            d.capping_required = "No"
            d.wall_to_internal = "No"
            d.invoke_share_topology = "Yes"
            d.non_conformal = "No"
            d()
        schritt("Geometrie beschreiben (Fluid + Solid)", beschreiben)
        names = wt.task_names()
        print("    Aufgaben:", names)

        # 4 --- Share Topology (macht die gemeinsamen Flaechen Luft/Bauteil konform)
        if "apply_share_topology" in names:
            schritt("Share Topology anwenden", lambda: wt.apply_share_topology())
        else:
            print("\n=== Share Topology: Aufgabe nicht vorhanden, uebersprungen")

        # 5 --- Raender
        schritt("Raender aktualisieren", lambda: task(wt, "update_boundaries")())

        # 6 --- Regionen typisieren: Luftraum = fluid, Bauteile = solid.
        # Bei Setup-Typ "fluid and solid" existieren die 5 Regionen schon (eine je Koerper);
        # create_regions (Fluid-Extraktion) wird NICHT benoetigt und schlaegt dort fehl.
        # Die aktuellen Regionen stehen in region_current_list/region_current_type_list.
        def regionen():
            ur = task(wt, "update_regions")
            # Setzen ueber region_type_list wirkt in 2024 R2 nicht -> Regionstyp per TUI am
            # Mesh-Objekt aendern. WICHTIG: ur.arguments() erst NACH der TUI-Aenderung lesen,
            # sonst haelt der Task den alten Zustand (alle solid) und schreibt ihn zurueck.
            fluid = FLUID
            meshing.tui.objects.volumetric_regions.change_type(
                STEP.stem, "(" + " ".join(fluid) + ")", "fluid")
            print("    Regionen (gesetzt per TUI):", fluid, "-> fluid")
            ur()
            danach = ur.arguments()
            n2 = list(danach.get("region_current_list") or [])
            t2 = list(danach.get("region_current_type_list") or [])
            print("    Regionen (danach):", list(zip(n2, t2)))
            if sorted(n for n, t in zip(n2, t2) if t == "fluid") != FLUID:
                raise RuntimeError("Nach Update Regions ist luftraum-luftraum NICHT die einzige Fluid-Region")
            return [n for n, t in zip(n2, t2) if t == "fluid"]
        fluid_regionen = schritt("Regionen festlegen", regionen)

        # 7a --- Prismenschichten (Grenzschicht) an den Fluid-Waenden. Bei Fehler: weiter ohne Prismen.
        def grenzschicht():
            bl = task(wt, "add_boundary_layer", "add_boundary_layers")
            try:
                print("    BL-Argumente:", list(bl.arguments().keys()))
            except Exception:
                pass
            # Prismen im FLUID (nicht in den Festkoerpern), nur an Waenden; per set_state
            bl.arguments.set_state({
                "add_child": "yes",
                "offset_method_type": "smooth-transition",
                "number_of_layers": BL_SCHICHTEN,
                "transition_ratio": BL_UEBERGANG,
                "rate": BL_WACHSTUM,
                "face_scope": {"regions_type": "fluid-regions", "grow_on": "only-walls"},
                "region_scope": fluid_regionen,
            })
            a = bl.arguments()
            print("    BL-Parameter:", {k: a.get(k) for k in ("offset_method_type", "number_of_layers",
                  "transition_ratio", "rate", "face_scope", "region_scope")})
            if a.get("number_of_layers") != BL_SCHICHTEN:
                raise RuntimeError("number_of_layers wurde nicht uebernommen")
            kind = bl.add_child_and_update(defer_update=False)
            kinder = []
            try:
                kinder = [t.name() for t in bl.tasks()]
            except Exception:
                pass
            print("    BL-Kind-Aufgaben:", kinder, "| letztes Kind:", kind)
            if not kinder and kind is None:
                raise RuntimeError("Kein BL-Kind angelegt - Grenzschicht NICHT definiert")
            return kinder
        bl_status = {"ok": False}
        def grenzschicht_wrapper():
            bl_status["kinder"] = grenzschicht()
            bl_status["ok"] = True
        try:
            schritt("Prismenschichten (Grenzschicht)", grenzschicht_wrapper)
        except Exception:
            print("    >>> Prismen fehlgeschlagen - Netz wird OHNE Prismenschichten erzeugt.", flush=True)

        # 7b --- Volumennetz
        def volumen():
            v = task(wt, "create_volume_mesh", "create_volume_mesh_wtm", "generate_the_volume_mesh_wtm")
            try:
                print("    Erlaubte volume_fill-Werte:", v.volume_fill.allowed_values())
            except Exception:
                pass
            v.arguments.set_state({"volume_fill": "poly-hexcore", "volume_fill_controls": {
                "hex_max_cell_length": s["hex_max"], "growth_rate": WACHSTUM}})
            a = v.arguments()
            print("    Volumennetz-Parameter:", a.get("volume_fill"),
                  {k: a.get("volume_fill_controls", {}).get(k) for k in ("hex_max_cell_length", "growth_rate")})
            v()
        schritt("Volumennetz erzeugen", volumen)

        # 7c --- Qualitaet verbessern (falls Aufgabe vorhanden)
        def verbessern():
            # Workflow-Task fehlt in 2024 R2 -> Root-Befehl ImproveVolumeMesh des Meshing-Datenmodells
            cmd = None
            for pfad in (lambda: meshing.meshing.ImproveVolumeMesh,
                         lambda: wt.improve_volume_mesh):
                try:
                    cmd = pfad()
                    break
                except Exception:
                    continue
            if cmd is None:
                print("    ImproveVolumeMesh nicht verfuegbar, uebersprungen")
                return
            vor = {"VIQualityIterations": 20, "VIQualityMinAngle": 0, "VIgnoreFeature": "yes"}
            # abwechselnd Skewness (Gate <= 0.85) und Orthogonal Quality (Gate > 0.15)
            for runde in range(3):
                for methode, grenze in (("Skewness", SKEWNESS_GRENZE - 0.01), ("Orthogonal", OQ_GRENZE + 0.01)):
                    cmd(QualityMethod=methode, CellQualityLimit=grenze, VMImprovePreferences=vor)
                    print(f"    Verbesserung Runde {runde + 1}: {methode} {grenze}")
        try:
            schritt("Volumennetz verbessern", verbessern)
        except Exception:
            print("    >>> Verbesserung fehlgeschlagen - Netz bleibt unveraendert.", flush=True)

        # 7d --- Zelltyp absichern: Luftraum = fluid (Workflow-Zustand kann solid bleiben)
        def zelltyp():
            for z in FLUID:
                meshing.tui.mesh.manage.type(z, "fluid")
            meshing.tui.mesh.manage.list()   # Tabelle mit Zellarten (wedge/poly-prism = Prismen)
        schritt("Zelltyp Luftraum = fluid", zelltyp)

        # 8 --- Speichern und Qualitaet
        # write_mesh ueberschreibt eine vorhandene Datei im Batch NICHT (Rueckfrage wird verneint)
        # -> alte Datei vorher entfernen und danach pruefen, dass sie neu geschrieben wurde.
        def speichern():
            t_start = time.time()
            if mesh_datei.exists():
                mesh_datei.unlink()
            meshing.tui.file.write_mesh(mesh_datei.as_posix())
            if not mesh_datei.exists() or mesh_datei.stat().st_mtime < t_start - 1:
                raise RuntimeError(f"Netzdatei wurde nicht geschrieben: {mesh_datei}")
            print(f"    geschrieben: {mesh_datei} ({mesh_datei.stat().st_size} Byte)")
        schritt(f"Netz speichern: {mesh_datei.name}", speichern)
        try:
            meshing.tui.mesh.check_quality()
        except Exception as e:
            print("    Hinweis: Qualitaetspruefung nicht moeglich:", e)
        try:
            print("\n=== Endkontrolle Skewness (alle Zellzonen)", flush=True)
            meshing.tui.report.quality_method("skewness")
            meshing.tui.report.cell_quality_limits("*", "()")
        except Exception as e:
            print("    Hinweis: Skewness-Endkontrolle nicht moeglich:", str(e)[:200])

        print("\n>>> PRISMENSCHICHTEN IM NETZ:", "JA " + str(bl_status.get("kinder")) if bl_status["ok"] else "NEIN (nicht definiert / fehlgeschlagen)")
        print(f"\nFertig. Netz: {mesh_datei}")
    finally:
        try:
            meshing.exit()
        except Exception:
            pass


if __name__ == "__main__":
    main()
