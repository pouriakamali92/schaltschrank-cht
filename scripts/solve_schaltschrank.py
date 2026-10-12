"""
solve_schaltschrank.py - CHT-Loesung des Schaltschranks mit Fluent (PyFluent 0.42.2):
natuerliche Konvektion innen + S2S-Strahlung innen + Waermeleitung in den Bauteilen,
Aussenwand als 2-mm-Stahlschale mit gemischter Randbedingung (Konvektion h + Strahlung).

Aufruf:
  python solve_schaltschrank.py [--netz grob|mittel|fein] [--fall NAME] [--nproc N]
                                [--set schluessel=wert ...]
Schluessel fuer --set (Standard aus config.py):
  eps_innen, eps_bauteil, eps_aussen, kfaktor (Faktor auf alle Bauteil-k), T_umg (degC),
  hfaktor (Faktor auf Aussen-h), h_fest (fester Aussen-h, schaltet Nachfuehrung ab),
  turb (laminar | k-omega-sst), dichte (incompressible-ideal-gas | boussinesq),
  P_Umrichter/P_Netzteil/P_Schuetz/P_Trafo (W), iter (max. Iterationen), fpsc (S2S-Cluster)
Ausgaben:
  laeufe/<fall>/  : Case/Data, ergebnis.json, monitor.csv, ebene_*.csv (Postprocessing-Daten)
  logs/solve_<fall>.txt : kurze Zusammenfassung (Zahlen, Pruefungen, Fehler)
"""
import argparse
import json
import math
import sys
import time
import traceback
from pathlib import Path

import ansys.fluent.core as pyfluent

import config as C
import handrechnung as HR

HIER = Path(__file__).resolve().parent
LAEUFE = HIER / "laeufe"
LOGS = HIER / "logs"

BAUTEILE = list(C.KOMPONENTEN)                      # Umrichter, Netzteil, Schuetz, Trafo
ZONE = {n: f"{n.lower()}-{n.lower()}" for n in BAUTEILE}
LUFT = "luftraum-luftraum"
AUSSEN_ROH = "luftraum-luftraum:1"
SEITEN = ["boden", "decke", "links", "rechts", "rueck", "tuer"]
# S2S-Clustering: Flaechen je Cluster so, dass die physikalische Clustergroesse ueber die
# Netzstufen etwa gleich bleibt (Verfeinerungsfaktor ~1.4 je Stufe -> Flaechen x2)
FPSC = {"grob": 2, "mittel": 4, "fein": 8}


def K(c):
    return c + 273.15


# ---------------------------------------------------------------------------
# Parameter
# ---------------------------------------------------------------------------
def standard_parameter(netz="mittel"):
    p = dict(
        netz=netz,
        eps_innen=C.EPS_WAND_INNEN, eps_bauteil=C.EPS_BAUTEIL, eps_aussen=C.EMISSIVITAET,
        kfaktor=1.0, T_umg=C.T_UMGEBUNG, hfaktor=1.0, h_fest=None,
        turb=C.TURBULENZ, dichte=C.DICHTE, iter=C.ITER_MAX, fpsc=FPSC.get(netz, 4),
        layout="A", luefter=0, rad_iter=50,
    )
    for n in BAUTEILE:
        p[f"P_{n}"] = C.KOMPONENTEN[n]
    return p


def parse_set(p, liste):
    for eintrag in liste or []:
        k, v = eintrag.split("=", 1)
        if k not in p:
            raise SystemExit(f"Unbekannter Parameter: {k}")
        if v.lower() in ("none", ""):
            p[k] = None
            continue
        try:
            p[k] = int(v) if k in ("iter", "fpsc") else float(v)
        except ValueError:
            p[k] = v
    return p


# ---------------------------------------------------------------------------
# Hilfen
# ---------------------------------------------------------------------------
class Protokoll:
    def __init__(self, datei):
        self.datei = datei
        self.zeilen = []

    def __call__(self, text=""):
        print(text, flush=True)
        self.zeilen.append(str(text))
        self.datei.write_text("\n".join(self.zeilen) + "\n", encoding="utf-8")


def berechne(st, namen):
    """Report-Definitionen auswerten -> {name: wert oder dict(zone: wert)}."""
    roh = st.solution.report_definitions.compute(report_defs=list(namen))
    out = {}
    for eintrag in roh:
        for name, wert in eintrag.items():
            if "(" in name:
                continue
            # Flux-Reports: Hauptwert ist "Net" = Fluss + Quellen im Gebiet!
            # Der eigentliche Waermestrom steht unter "<name>(without-sources)".
            wert = eintrag.get(f"{name}(without-sources)", wert)
            out[name] = wert[0] if isinstance(wert, (list, tuple)) else wert
            if f"{name}(User Energy Source)" in eintrag:
                out["_quelle_fluent"] = eintrag[f"{name}(User Energy Source)"][0]
    return out


def vol_report(st, name, typ, feld, zonen):
    rd = st.solution.report_definitions.volume
    rd.create(name=name)
    o = rd[name]
    o.report_type = typ
    if feld:
        o.field = feld
    o.cell_zones = zonen
    return name


def flaechen_report(st, name, typ, feld, flaechen):
    rd = st.solution.report_definitions.surface
    rd.create(name=name)
    o = rd[name]
    o.report_type = typ
    if feld:
        o.field = feld
    o.surface_names = flaechen
    return name


def fluss_report(st, name, raender):
    rd = st.solution.report_definitions.flux
    rd.create(name=name)
    o = rd[name]
    o.report_type = "flux-heattransfer"
    o.boundaries = raender
    return name


# ---------------------------------------------------------------------------
# Setup
# ---------------------------------------------------------------------------
def starte(nproc, cwd=None):
    """Fluent-Solver starten. Lizenz: max. 12 Kerne je Lauf (16 -> 'license not available')."""
    return pyfluent.launch_fluent(mode="solver", dimension=3, precision="double",
                                  processor_count=nproc, ui_mode="no_gui", cleanup_on_exit=True,
                                  cwd=str(cwd) if cwd else None)


def netz_lesen(s, netzdatei, log, solid_soll=None):
    st = s.settings
    st.file.read_mesh(file_name=Path(netzdatei).as_posix())
    czc = st.setup.cell_zone_conditions
    fluid = sorted(czc.fluid.get_object_names())
    if len(fluid) > 1 and all(z == LUFT or z.startswith("oeffnung_") for z in fluid):
        # Luefter-Variante: Oeffnungskaesten mit dem Luftraum zu EINER Fluid-Zone verschmelzen
        st.mesh.modify_zones.merge_zones(zone_names=fluid)
        fluid_neu = sorted(czc.fluid.get_object_names())
        log(f"Fluid-Zonen {fluid} verschmolzen -> {fluid_neu}")
        if len(fluid_neu) == 1 and fluid_neu[0] != LUFT:
            st.mesh.modify_zones.zone_name(zone_name=fluid_neu[0], new_name=LUFT)
        fluid = sorted(czc.fluid.get_object_names())
        # Die frueheren Grenzflaechen Luftraum/Oeffnungskasten sind Wand+Schatten-Paare -> durchlaessig machen
        W = st.setup.boundary_conditions.wall
        innen = [w for w in W.get_object_names() if w.startswith(f"{LUFT}-oeffnung_") and not w.endswith("-shadow")]
        for w in innen:
            st.setup.boundary_conditions.set_zone_type(zone_list=[w], new_type="interior")
        rest = [w for w in W.get_object_names() if w.startswith(f"{LUFT}-oeffnung_")]
        log(f"Innenwaende Luftraum/Oeffnungskasten -> interior: {innen}; verbleibend als Wand: {rest}")
        if rest:
            raise RuntimeError(f"Grenzflaechen zu den Oeffnungskaesten noch Waende: {rest}")
    solid = sorted(czc.solid.get_object_names())
    log(f"Zonen: fluid={fluid} solid={solid}")
    if fluid != [LUFT] or set(solid) != set(solid_soll or ZONE.values()):
        raise RuntimeError("Zonenpruefung fehlgeschlagen (Luftraum muss einzige Fluid-Zone sein)")
    # Aussenwand (eine Zone) nach Flaechenwinkel in die 6 Schrankseiten trennen
    s.tui.mesh.modify_zones.sep_face_zone_angle(AUSSEN_ROH, 60, "yes")
    waende = st.setup.boundary_conditions.wall.get_object_names()
    aussen = sorted(w for w in waende if w == AUSSEN_ROH or w.startswith(AUSSEN_ROH + ":"))
    log(f"Aussenwand-Zonen nach Trennung: {aussen}")
    return aussen


def modelle_und_stoffe(s, p, log):
    st = s.settings
    og = st.setup.general.operating_conditions
    og.gravity.enable = True
    og.gravity.components = [0.0, -C.G, 0.0]
    og.operating_pressure = 101325.0

    st.setup.models.energy.enabled = True
    vis = st.setup.models.viscous
    if p["turb"] == "laminar":
        vis.model = "laminar"
    else:
        vis.model = "k-omega"
        vis.k_omega_model = "sst"
    # S2S erst in s2s_sichtfaktoren() aktivieren (Initialisierung braucht sonst eine VF-Datei)

    luft = st.setup.materials.fluid["air"]
    T0 = K(p["T_umg"] + 10.0)           # Bezugstemperatur Boussinesq (ca. Innenluft)
    if p["dichte"] == "boussinesq":
        luft.density.option = "boussinesq"
        luft.density.value = 101325.0 / (287.05 * T0)
        luft.therm_exp_coeff.value = 1.0 / T0
        og.operating_temperature = T0
        log(f"Boussinesq: rho0 = {101325.0 / (287.05 * T0):.4f} kg/m3, beta = 1/{T0:.1f} K, T0 = {T0:.1f} K")
    else:
        luft.density.option = "incompressible-ideal-gas"
    luft.viscosity.option = "sutherland"
    luft.thermal_conductivity.value = C.LUFT_LAMBDA
    luft.specific_heat.value = 1006.43

    for n, m in C.BAUTEIL_MATERIAL.items():
        name = f"ersatz-{n.lower()}"
        st.setup.materials.solid.create(name)
        mat = st.setup.materials.solid[name]
        mat.density.value = m["rho"]
        mat.specific_heat.value = m["cp"]
        mat.thermal_conductivity.value = m["k"] * p["kfaktor"]
    st.setup.materials.solid.create("stahl")
    stahl = st.setup.materials.solid["stahl"]
    stahl.density.value = C.WAND_RHO
    stahl.specific_heat.value = C.WAND_CP
    stahl.thermal_conductivity.value = C.WAND_K


def bauteile_setzen(s, p, log, nur=None):
    st = s.settings
    czc = st.setup.cell_zone_conditions
    for n in (nur or BAUTEILE):
        dx, dy, dz = C.BLOCK_SIZE[n]
        V = dx * dy * dz
        q = p[f"P_{n}"] / V
        z = czc.solid[ZONE[n]]
        z.general.material = f"ersatz-{n.lower()}"
        z.sources.enable = True
        z.sources.terms = {"energy": [{"option": "value", "value": q}]}
        log(f"  {n:10s}: P = {p[f'P_{n}']:6.1f} W, V = {V * 1e3:.3f} dm3, q''' = {q:9.1f} W/m3, "
            f"k = {C.BAUTEIL_MATERIAL[n]['k'] * p['kfaktor']:.2f} W/mK")


def seiten_zuordnen(s, aussen):
    """Schrankseite je Aussenwand-Zone aus dem Flaechenschwerpunkt (benoetigt Initialisierung)."""
    st = s.settings
    namen = []
    for w in aussen:
        for a in "xyz":
            namen.append(flaechen_report(st, f"c{a}-{w.replace(':', '_')}", "surface-areaavg",
                                         f"{a}-coordinate", [w]))
    werte = berechne(st, namen)
    zuordnung = {}
    eps = 0.01
    for w in aussen:
        x, y, z = (werte[f"c{a}-{w.replace(':', '_')}"] for a in "xyz")
        if y < eps:
            seite = "boden"
        elif y > C.H - eps:
            seite = "decke"
        elif x < eps:
            seite = "links"
        elif x > C.B - eps:
            seite = "rechts"
        elif z < eps:
            seite = "rueck"
        elif z > C.T - eps:
            seite = "tuer"
        else:
            raise RuntimeError(f"Aussenwand {w} keiner Seite zuordenbar: ({x:.3f}, {y:.3f}, {z:.3f})")
        zuordnung.setdefault(seite, []).append(w)
    for n in namen:
        del st.solution.report_definitions.surface[n]
    fehlend = [x for x in SEITEN if x not in zuordnung]
    if fehlend:
        raise RuntimeError(f"Schrankseiten ohne Zone: {fehlend}")
    return zuordnung


def luefter_einbauen(s, p, seiten, log):
    """Erzwungene Kuehlung. Geometrie: build_geometry.py --luefter (Oeffnungskaesten in der Tuer, im Solver
    mit dem Luftraum verschmolzen; ihre Aussenflaeche ist eine eigene Randzone "oeffnung_*").
    (Ein Heraustrennen per Zellregister + sep_face_zone_mark markierte in 2024 R2 keine Zellen -> verworfen.)
    Eintritt = intake-fan (Luefterkennlinie minus Eintrittsfilter als Drucksprung-Polynom in der
    Normalgeschwindigkeit), Austritt = outlet-vent (Filter als Verlustbeiwert).
    Rueckgabe: {"ein": zone, "aus": zone, "A_ein": m2, "A_aus": m2, ...}."""
    st = s.settings
    W = st.setup.boundary_conditions.wall
    out = {}
    for key in ("ein", "aus"):
        z = [w for w in W.get_object_names() if w.startswith(f"oeffnung_{key}") and not w.endswith("-shadow")]
        if len(z) != 1:
            raise RuntimeError(f"Oeffnung {key}: erwartet genau 1 Randzone, gefunden {z} "
                               f"(Netz aus build_geometry.py --luefter?)")
        out[key] = z[0]
    # Flaechen der tatsaechlich herausgetrennten Facetten (netzabhaengig) bestimmen
    namen = [flaechen_report(st, f"a-{k}", "surface-area", None, [out[k]]) for k in ("ein", "aus")]
    werte = berechne(st, namen)
    out["A_ein"], out["A_aus"] = werte["a-ein"], werte["a-aus"]
    for n in namen:
        del st.solution.report_definitions.surface[n]
    bc = st.setup.boundary_conditions
    bc.set_zone_type(zone_list=[out["ein"]], new_type="intake-fan")
    bc.set_zone_type(zone_list=[out["aus"]], new_type="outlet-vent")
    A_e = out["A_ein"]
    rho = 101325.0 / (287.05 * K(p["T_umg"]))
    # dp(v) = DP0 (1 - v A / V0) - zeta rho/2 v^2   (Luefterkennlinie minus Eintrittsfilter)
    koeff = [C.LUEFTER_DP0, -C.LUEFTER_DP0 * A_e / C.LUEFTER_V0, -C.FILTER_ZETA * rho / 2]
    ein = bc.intake_fan[out["ein"]]
    ein.momentum.pressure_jump = {"option": "polynomial", "function_of": "normal-velocity",
                                  "polynomial": {"function_of": "normal-velocity", "coefficients": koeff}}
    ein.thermal.total_temperature.value = K(p["T_umg"])
    aus = bc.outlet_vent[out["aus"]]
    aus.momentum.loss_coefficient = {"option": "polynomial", "function_of": "normal-velocity",
                                     "polynomial": {"function_of": "normal-velocity",
                                                    "coefficients": [C.FILTER_ZETA]}}
    aus.thermal.backflow_total_temperature.value = K(p["T_umg"])
    log(f"Luefter: Eintritt {out['ein']} ({1e4 * A_e:.1f} cm2), Austritt {out['aus']} "
        f"({1e4 * out['A_aus']:.1f} cm2), Drucksprung-Polynom {[round(k, 3) for k in koeff]}")
    out["koeff"] = koeff
    out["rho"] = rho
    return out


def luefter_strahlung(s, p, oeffn, log):
    """Pruefung: Oeffnungen strahlen als schwarzer Strahler mit der Randtemperatur (= T_umg, Fluent-Vorgabe
    'Boundary Temperature'); nur protokollieren."""
    bc = s.settings.setup.boundary_conditions
    for typ, z in (("intake_fan", oeffn["ein"]), ("outlet_vent", oeffn["aus"])):
        rad = getattr(bc, typ)[z].radiation()
        log(f"  Strahlung Oeffnung {z}: {rad.get('external_black_body_temperature_method')}, "
            f"eps {rad.get('internal_emissivity', {}).get('value')}")


def luefter_betriebspunkt(oeffn):
    """Analytischer Betriebspunkt: Luefter-Kennlinie mit Eintrittsfilter = Austrittsfilter
    (Innenverluste und Kamineffekt vernachlaessigt) -> Vergleich mit CFD."""
    c0, c1, c2 = oeffn["koeff"]
    A_e, A_a, rho = oeffn["A_ein"], oeffn["A_aus"], oeffn["rho"]
    # c0 + c1 v + c2 v^2 = zeta rho/2 (v A_e/A_a)^2  (gleicher Volumenstrom, v = Eintrittsgeschw.)
    a2 = c2 - C.FILTER_ZETA * rho / 2 * (A_e / A_a) ** 2
    v = (-c1 - (c1 * c1 - 4 * a2 * c0) ** 0.5) / (2 * a2)
    return v * A_e, c0 + c1 * v + c2 * v * v


def h_start(p):
    """Aussen-h (konvektiv) je Seite aus der Handrechnung (gleichmaessige Wandtemperatur)."""
    if p["h_fest"] is not None:
        return {x: p["h_fest"] for x in SEITEN if x != "boden"}
    Tw = HR.loese_wandtemperatur(sum(p[f"P_{n}"] for n in BAUTEILE), p["T_umg"], p["eps_aussen"], p["hfaktor"])
    h = {}
    for x in SEITEN:
        if x == "boden":
            continue
        h[x] = p["hfaktor"] * HR.h_konvektion(Tw, "deckel" if x == "decke" else "senkrecht", p["T_umg"])
    return h


def aussenwand_setzen(s, p, seiten, h):
    st = s.settings
    W = st.setup.boundary_conditions.wall
    for seite, zonen in seiten.items():
        for w in zonen:
            th = W[w].thermal
            if seite == "boden" and C.BODEN_ADIABAT:
                th.thermal_condition = "Heat Flux"
                th.heat_flux.value = 0.0
            else:
                th.thermal_condition = "Mixed"
                th.heat_transfer_coeff.value = h[seite]
                th.free_stream_temp.value = K(p["T_umg"])
                th.external_emissivity.value = p["eps_aussen"]
                th.ext_rad_temperature.value = K(p["T_umg"])
            th.material = "stahl"
            th.wall_thickness.value = C.WAND
            th.enable_shell_conduction = True
            W[w].radiation.internal_emissivity.value = p["eps_innen"]


def innenwaende_setzen(s, p, aussen, log):
    """Emissivitaet der Bauteiloberflaechen (luftseitige Wand); erst nach Aktivieren von S2S."""
    st = s.settings
    W = st.setup.boundary_conditions.wall
    gesetzt, fehler = 0, []
    for w in W.get_object_names():
        # nur Grenzflaechen Luft/Bauteil (Namen "luftraum-luftraum-<bauteil>..."); Festkoerper-
        # Festkoerper-Grenzflaechen (z.B. Kern/Wicklung) strahlen nicht
        if w in aussen or "luftraum" not in w:
            continue
        try:
            W[w].radiation.internal_emissivity.value = p["eps_bauteil"]
            gesetzt += 1
        except Exception as e:
            fehler.append(w)   # Schattenseite im Festkoerper: keine Strahlungsparameter aktiv
    luftseitig_fehlt = [w for w in fehler if not w.endswith("-shadow")]
    log(f"Emissivitaet Bauteile {p['eps_bauteil']} auf {gesetzt} Waenden gesetzt"
        + (f"; NICHT gesetzt (luftseitig!): {luftseitig_fehlt}" if luftseitig_fehlt else ""))
    if luftseitig_fehlt:
        raise RuntimeError("Emissivitaet auf luftseitigen Bauteilwaenden nicht setzbar")


def numerik(s, p):
    st = s.settings
    m = st.solution.methods
    m.p_v_coupling.flow_scheme = "Coupled"
    m.discretization_scheme["pressure"] = "body-force-weighted"
    for gl in ("mom", "temperature", "k", "omega"):
        try:
            m.discretization_scheme[gl] = "second-order-upwind"
        except Exception:
            pass
    # Abbruch steuert das Skript ueber Monitorgroessen -> Residuenkriterien sehr klein
    eqs = st.solution.monitor.residual.equations
    for gl in eqs.get_object_names() if hasattr(eqs, "get_object_names") else eqs().keys():
        try:
            eqs[gl].absolute_criteria = 1e-12
        except Exception:
            pass


def s2s_sichtfaktoren(s, p, ordner, log):
    st = s.settings
    s2s = st.setup.models.radiation.s2s
    s2s.clustering_settings.faces_per_cluster = {"option": "manual",
                                                 "global_faces_per_surface_cluster": int(p["fpsc"])}
    if p.get("rad_iter"):
        # Radiositaetsgleichungen je Fluss-Iteration genau loesen. Fluent-Vorgabe (5 It., 1e-3) ergab bei
        # eps 0,6 einen bleibenden Bilanzfehler von 1,9 %; mit 50 It. / 1e-6: 0,13 % (docs/log.md)
        s2s.radiosity_solver_control.maximum_radiation_iterations = int(p["rad_iter"])
        s2s.radiosity_solver_control.residual_convergence_criteria = 1e-6
    log(f"S2S: Sichtfaktoren {s2s.viewfactor_settings()}, Radiositaet {s2s.radiosity_solver_control()}")
    vf = (ordner / "sichtfaktoren.s2s.h5").as_posix()
    t0 = time.time()
    try:
        s2s.compute_write_vf(file_name=vf)
    except Exception as e:
        log(f"  compute_write_vf (settings) fehlgeschlagen: {str(e)[:150]} -> TUI")
        s.tui.define.models.radiation.s2s_parameters.compute_write_vf(vf)
    log(f"S2S-Sichtfaktoren berechnet ({time.time() - t0:.0f} s, {p['fpsc']} Flaechen/Cluster)")


# ---------------------------------------------------------------------------
# Loesen
# ---------------------------------------------------------------------------
def reports_anlegen(s, seiten, aussen, gruppen=None):
    """gruppen: Bezeichnung -> Liste von Zellzonen (Standard: je Bauteil eine Zone)."""
    st = s.settings
    r = {"tmax": {}, "tavg": {}, "q_bt": {}, "t_seite": {}, "q_seite": {}}
    W = st.setup.boundary_conditions.wall.get_object_names()
    for n, zonen in (gruppen or {n: [ZONE[n]] for n in BAUTEILE}).items():
        r["tmax"][n] = vol_report(st, f"tmax-{n.lower()}", "volume-max", "temperature", zonen)
        r["tavg"][n] = vol_report(st, f"tavg-{n.lower()}", "volume-average", "temperature", zonen)
        # Waermestrom Bauteil -> Luft ueber die luftseitigen Schattenwaende
        schatten = [w for w in W if w.endswith("-shadow") and w.startswith(LUFT)
                    and any(f"-{z}" in w for z in zonen)]
        r["q_bt"][n] = fluss_report(st, f"q-{n.lower()}", schatten)
    r["tavg_luft"] = vol_report(st, "tavg-luft", "volume-average", "temperature", [LUFT])
    r["tmax_luft"] = vol_report(st, "tmax-luft", "volume-max", "temperature", [LUFT])
    r["vmax_luft"] = vol_report(st, "vmax-luft", "volume-max", "velocity-magnitude", [LUFT])
    r["q_aussen"] = fluss_report(st, "q-aussen", aussen)
    for x, zonen in seiten.items():
        r["t_seite"][x] = flaechen_report(st, f"t-{x}", "surface-areaavg", "temperature", zonen)
        r["q_seite"][x] = fluss_report(st, f"q-{x}", zonen)
    return r


def alle_namen(r):
    out = []
    for v in r.values():
        out.extend(v.values() if isinstance(v, dict) else [v])
    return out


def loese(s, p, r, seiten, h, log, ordner, Q=None, monitor="monitor.csv", chunk=50, fenster=1000,
          it_min=1000):
    """Pseudo-transient iterieren, bis die Monitorgroessen quasi-stationaer sind; Aussen-h nachfuehren.

    Natuerliche Konvektion bei Ra ~ 1e9 ist leicht instationaer (kleine Schwankungen der Fahnen).
    Kriterium daher ueber ein Fenster (fenster Iterationen, Abtastung alle chunk Iterationen):
      Drift = |Mittel 2. Haelfte - Mittel 1. Haelfte| < max(0.1 K, 0.5 x Schwankungsamplitude)
      fuer alle T-Monitore der Bauteile/Waende/Luftmittel (das lokale Luftmaximum schwankt stark und
      wird nur berichtet), Drift Q_aussen < 0.5 %, Energiebilanz des Fenstermittels < 1 %,
      mindestens it_min Iterationen. Die Schwankungen der Fahnen haben Perioden von ~1000 Iterationen;
      ein kuerzeres Fenster misst die Schwingung selbst statt einer Drift (docs/log.md).
    Ergebnis = Fenstermittelwerte; Schwankungsamplitude (max-min)/2 wird mitgeschrieben.
    """
    st = s.settings
    namen = alle_namen(r)
    tnamen = [n for n in namen if n.startswith(("tmax-", "tavg-", "t-")) and n != "tmax-luft"]
    labels = list(r["tmax"])
    Q = Q if Q is not None else sum(p[f"P_{n}"] for n in BAUTEILE)
    hist = []
    it = 0
    fenster_start = 0           # Iteration, ab der das aktuelle Fenster zaehlt (nach h-Update neu)
    h_iter = 0
    nachfuehren = p["h_fest"] is None
    konvergiert = False
    mittel, ampl = None, None
    t0 = time.time()
    n_fenster = fenster // chunk
    while it < p["iter"]:
        st.solution.run_calculation.iterate(iter_count=chunk)
        it += chunk
        v = berechne(st, namen)
        v["iter"] = it
        hist.append(v)
        if it % 100 == 0:
            bil = (Q + v["q-aussen"]) / Q
            tm = {n: v[r["tmax"][n]] - 273.15 for n in labels}
            log(f"  it {it:5d}: Tmax " + " ".join(f"{n[:4]}={tm[n]:6.2f}" for n in labels)
                + f" | Tluft={v['tavg-luft'] - 273.15:6.2f} | Q_aussen={v['q-aussen']:8.2f} W"
                + f" | Bilanz={100 * bil:+6.2f} % | {time.time() - t0:5.0f} s")
        w = [x for x in hist if x["iter"] > fenster_start][-n_fenster:]
        if len(w) < n_fenster or w[0]["iter"] <= fenster_start:
            continue
        hz = len(w) // 2
        m1 = {n: sum(x[n] for x in w[:hz]) / hz for n in namen}
        m2 = {n: sum(x[n] for x in w[hz:]) / (len(w) - hz) for n in namen}
        mittel = {n: sum(x[n] for x in w) / len(w) for n in namen}
        ampl = {n: (max(x[n] for x in w) - min(x[n] for x in w)) / 2 for n in namen}
        drift_T = max(abs(m2[n] - m1[n]) for n in tnamen)
        drift_ok = all(abs(m2[n] - m1[n]) < max(0.1, 0.5 * ampl[n]) for n in tnamen)
        drift_Q = abs(m2["q-aussen"] - m1["q-aussen"]) / Q
        bil = (Q + mittel["q-aussen"]) / Q
        if not (drift_ok and drift_Q < 0.005 and abs(bil) < 0.01 and it - fenster_start >= it_min):
            continue
        # Aussen-h mit CFD-Wandtemperaturen (Fenstermittel) nachfuehren
        if nachfuehren and h_iter < C.H_ITER_MAX:
            h_neu = {}
            for x in h:
                Tw = mittel[r["t_seite"][x]] - 273.15
                lage = "deckel" if x == "decke" else "senkrecht"
                h_neu[x] = p["hfaktor"] * HR.h_konvektion(Tw, lage, p["T_umg"])
            aend = max(abs(h_neu[x] - h[x]) / h[x] for x in h)
            log(f"  h-Nachfuehrung {h_iter + 1} (it {it}): max. Aenderung {100 * aend:.1f} %  "
                + " ".join(f"{x}={h_neu[x]:.2f}" for x in h_neu))
            if aend > C.H_ITER_TOL:
                h.update(h_neu)
                aussenwand_setzen(s, p, seiten, h)
                h_iter += 1
                fenster_start = it
                continue
        konvergiert = True
        log(f"Quasi-stationaer nach {it} Iterationen (Fenster {fenster}): Drift T max {drift_T:.4f} K, "
            f"Drift Q {100 * drift_Q:.3f} %, Bilanz Fenstermittel {100 * bil:+.3f} %")
        log("  Schwankungsamplituden (K bzw. W): " + ", ".join(
            f"{n}={ampl[n]:.3f}" for n in [r["tmax"][b] for b in labels] + ["tavg-luft", "q-aussen"]))
        break
    if mittel is None:          # kein volles Fenster: letzte Werte
        mittel, ampl = dict(hist[-1]), {n: float("nan") for n in namen}
    # Verlauf speichern
    with open(ordner / monitor, "w", encoding="utf-8") as f:
        f.write("iter," + ",".join(namen) + "\n")
        for x in hist:
            f.write(f"{x['iter']}," + ",".join(f"{x.get(n, float('nan')):.6g}" for n in namen) + "\n")
    return konvergiert, it, mittel, ampl, h


# ---------------------------------------------------------------------------
# Auswertung
# ---------------------------------------------------------------------------
def residuen(s):
    try:
        idx, daten = s.monitors.get_monitor_set_data(monitor_set_name="residual")
        return {k: float(v[-1]) for k, v in daten.items() if len(v)}
    except Exception as e:
        return {"fehler": str(e)[:100]}


def grenzschicht_pruefen(s, p, aussen, log):
    """Wandnaechste Zellhoehe (2 x Zentrumsabstand) und y+ auf den Waenden der Luft."""
    st = s.settings
    W = st.setup.boundary_conditions.wall.get_object_names()
    luftseite = [w for w in W if not w.endswith("-shadow")]
    bt_waende = [w for w in luftseite if w not in aussen]
    namen = [flaechen_report(st, "wd-aussen", "surface-areaavg", "cell-wall-distance", aussen),
             flaechen_report(st, "wd-bauteil", "surface-areaavg", "cell-wall-distance", bt_waende),
             flaechen_report(st, "wdmax-bauteil", "surface-facetmax", "cell-wall-distance", bt_waende)]
    if p["turb"] != "laminar":
        namen += [flaechen_report(st, "yp-avg", "surface-areaavg", "y-plus", luftseite),
                  flaechen_report(st, "yp-max", "surface-facetmax", "y-plus", luftseite)]
    try:
        v = berechne(st, namen)
    except Exception as e:
        log(f"  Grenzschichtpruefung nicht moeglich: {str(e)[:150]}")
        return {}
    out = {"h1_aussen_mm": 2e3 * v["wd-aussen"], "h1_bauteil_mm": 2e3 * v["wd-bauteil"],
           "h1_bauteil_max_mm": 2e3 * v["wdmax-bauteil"]}
    if "yp-avg" in v:
        out.update(yplus_mittel=v["yp-avg"], yplus_max=v["yp-max"])
    return out


def ebenen_exportieren(s, ordner, log):
    st = s.settings
    ps = st.results.surfaces.plane_surface
    ebenen = {
        "z060": ("xy-plane", "z", 0.06),   # schneidet alle Bauteile (Rueckwand + 0.04 m)
        "z300": ("xy-plane", "z", 0.30),   # Luftraum vor den Bauteilen
        "x200": ("yz-plane", "x", 0.20),   # durch den Umrichter
        "x650": ("yz-plane", "x", 0.65),   # durch den Trafo
    }
    for name, (methode, achse, wert) in ebenen.items():
        n = f"ebene-{name}"
        ps.create(name=n)
        ps[n].method = methode
        setattr(ps[n], achse, wert)
        datei = ordner / f"ebene_{name}.csv"
        try:
            st.file.export.ascii(file_name=datei.as_posix(), surface_name_list=[n], delimiter="comma",
                                 cell_func_domain=["temperature", "velocity-magnitude",
                                                   "x-velocity", "y-velocity", "z-velocity"],
                                 location="node")
        except Exception as e:
            log(f"  Export {n} fehlgeschlagen: {str(e)[:150]}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--netz", default="mittel")
    ap.add_argument("--fall", default=None)
    ap.add_argument("--nproc", type=int, default=C.NPROC_SOLVER)
    ap.add_argument("--set", nargs="*", default=[])
    a = ap.parse_args()
    p = parse_set(standard_parameter(a.netz), a.set)
    fall = a.fall or f"basis_{a.netz}"
    ordner = LAEUFE / fall
    ordner.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(exist_ok=True)
    log = Protokoll(LOGS / f"solve_{fall}.txt")
    log(f"Fall {fall}  ({time.strftime('%Y-%m-%d %H:%M')})")
    log("Parameter: " + json.dumps(p))
    lay = ("" if p["layout"] == "A" else f"_layout{p['layout']}") + ("_luefter" if p["luefter"] else "")
    netzdatei = HIER / "mesh" / f"schaltschrank{lay}_{a.netz}.msh.h5"
    t_start = time.time()
    s = starte(a.nproc, ordner)
    ergebnis = {"fall": fall, "parameter": p, "ok": False}
    try:
        aussen = netz_lesen(s, netzdatei, log)
        modelle_und_stoffe(s, p, log)
        log("Bauteile:")
        bauteile_setzen(s, p, log)
        numerik(s, p)
        # Initialisierung (fuer Schwerpunkte der Aussenflaechen), dann Aussen-RB
        st = s.settings
        st.solution.initialization.initialization_type = "standard"
        st.solution.initialization.defaults["temperature"] = K(p["T_umg"] + 10.0)
        st.solution.initialization.initialize()
        seiten = seiten_zuordnen(s, aussen)
        log("Seiten: " + json.dumps(seiten))
        oeffn = luefter_einbauen(s, p, seiten, log) if p["luefter"] else None
        h = h_start(p)
        log("Aussen-h Start (konvektiv, W/m2K): " + " ".join(f"{x}={v:.2f}" for x, v in h.items())
            + ("  | Boden adiabat" if C.BODEN_ADIABAT else ""))
        st.setup.models.radiation.model = "s2s"
        aussenwand_setzen(s, p, seiten, h)
        innenwaende_setzen(s, p, aussen, log)
        if oeffn:
            luefter_strahlung(s, p, oeffn, log)
        s2s_sichtfaktoren(s, p, ordner, log)
        raender = aussen + ([oeffn["ein"], oeffn["aus"]] if oeffn else [])
        r = reports_anlegen(s, seiten, raender)
        if oeffn:
            rd = st.solution.report_definitions.flux
            rd.create(name="m-ein")
            rd["m-ein"].report_type = "flux-massflow"
            rd["m-ein"].boundaries = [oeffn["ein"]]
            r["m_ein"] = "m-ein"
        st.solution.initialization.initialize()

        konv, it, v, ampl, h = loese(s, p, r, seiten, h, log, ordner)
        Q = sum(p[f"P_{n}"] for n in BAUTEILE)
        bil = (Q + v["q-aussen"]) / Q
        ergebnis.update(
            konvergiert=konv, iterationen=it, Q_quellen_W=Q, Q_aussen_W=-v["q-aussen"],
            Q_aussen_ampl_W=ampl["q-aussen"],
            Tmax_ampl_K={n: ampl[r["tmax"][n]] for n in BAUTEILE},
            Tluft_ampl_K=ampl["tavg-luft"],
            bilanzfehler_proz=100 * bil,
            Tmax_C={n: v[r["tmax"][n]] - 273.15 for n in BAUTEILE},
            Tmittel_C={n: v[r["tavg"][n]] - 273.15 for n in BAUTEILE},
            Q_bauteil_W={n: -v[r["q_bt"][n]] for n in BAUTEILE},
            Tluft_mittel_C=v["tavg-luft"] - 273.15, Tluft_max_C=v["tmax-luft"] - 273.15,
            vmax_luft=v["vmax-luft"],
            T_seite_C={x: v[r["t_seite"][x]] - 273.15 for x in SEITEN},
            Q_seite_W={x: -v[r["q_seite"][x]] for x in SEITEN},
            h_aussen=h, residuen=residuen(s),
        )
        if oeffn:
            V_cfd = abs(v["m-ein"]) / oeffn["rho"]
            V_an, dp_an = luefter_betriebspunkt(oeffn)
            ergebnis["luefter"] = dict(V_cfd_m3h=3600 * V_cfd, V_analytisch_m3h=3600 * V_an,
                                       dp_analytisch_Pa=dp_an, A_ein_m2=oeffn["A_ein"], A_aus_m2=oeffn["A_aus"],
                                       koeff=oeffn["koeff"])
            log(f"Luefter: Volumenstrom CFD {3600 * V_cfd:.1f} m3/h, analytischer Betriebspunkt "
                f"{3600 * V_an:.1f} m3/h bei {dp_an:.1f} Pa (ohne Innenverluste/Kamineffekt)")
        ergebnis["grenzschicht"] = grenzschicht_pruefen(s, p, aussen, log)
        ebenen_exportieren(s, ordner, log)
        st.file.write_case_data(file_name=(ordner / f"{fall}.cas.h5").as_posix())
        ergebnis["laufzeit_s"] = time.time() - t_start
        ergebnis["ok"] = True

        # ---- Zusammenfassung / Plausibilitaet
        log("")
        log("=== ERGEBNIS ===")
        log(f"Konvergiert (Monitore stationaer): {'JA' if konv else 'NEIN'} nach {it} Iterationen, "
            f"Laufzeit {ergebnis['laufzeit_s'] / 60:.1f} min")
        log(f"Energiebilanz: Quellen {Q:.2f} W, Abgabe aussen {-v['q-aussen']:.2f} W, "
            f"Fehler {100 * bil:+.3f} %  -> {'OK (<1 %)' if abs(bil) < 0.01 else 'VERLETZT'}")
        for n in BAUTEILE:
            log(f"  {n:10s}: Tmax {ergebnis['Tmax_C'][n]:6.2f} C, Tmittel {ergebnis['Tmittel_C'][n]:6.2f} C, "
                f"abgegeben {ergebnis['Q_bauteil_W'][n]:7.2f} W (Quelle {p[f'P_{n}']:.1f} W)")
        log(f"Luft: Tmittel {ergebnis['Tluft_mittel_C']:.2f} C, Tmax {ergebnis['Tluft_max_C']:.2f} C, "
            f"vmax {ergebnis['vmax_luft']:.3f} m/s")
        log("Wand (innen, flaechengemittelt) / Abgabe je Seite:")
        for x in SEITEN:
            log(f"  {x:7s}: T {ergebnis['T_seite_C'][x]:6.2f} C, Q {ergebnis['Q_seite_W'][x]:7.2f} W")
        Tw_hr = HR.loese_wandtemperatur(Q, p["T_umg"], p["eps_aussen"], p["hfaktor"])
        A_ab = sum(ergebnis["Q_seite_W"][x] for x in SEITEN if x != "boden")
        Tw_cfd = sum(ergebnis["T_seite_C"][x] * ergebnis["Q_seite_W"][x] for x in SEITEN if x != "boden") / A_ab
        ergebnis.update(Tw_hand_C=Tw_hr, Tw_cfd_C=Tw_cfd)
        log(f"Handrechnung: gleichmaessige Wandtemperatur {Tw_hr:.2f} C (Untergrenze fuer Innenluft); "
            f"CFD: Wand (Q-gewichtet) {Tw_cfd:.2f} C, Luft-Mittel {ergebnis['Tluft_mittel_C']:.2f} C "
            f"(+{ergebnis['Tluft_mittel_C'] - Tw_hr:.2f} K ueber Handrechnung)")
        log("Schwankung (Amplitude im Fenster): " + ", ".join(
            f"{n} {ergebnis['Tmax_ampl_K'][n]:.2f} K" for n in BAUTEILE)
            + f", Luft {ergebnis['Tluft_ampl_K']:.3f} K, Q_aussen {ergebnis['Q_aussen_ampl_W']:.2f} W")
        log(f"Grenzschicht: {json.dumps(ergebnis['grenzschicht'])}")
        log(f"Residuen (letzte): {json.dumps(ergebnis['residuen'])}")
    except Exception as e:
        log("!!! FEHLER: " + repr(e))
        log(traceback.format_exc())
        ergebnis["fehler"] = repr(e)
    finally:
        (ordner / "ergebnis.json").write_text(json.dumps(ergebnis, indent=2), encoding="utf-8")
        try:
            s.exit()
        except Exception:
            pass
    sys.exit(0 if ergebnis.get("ok") else 1)


if __name__ == "__main__":
    main()
