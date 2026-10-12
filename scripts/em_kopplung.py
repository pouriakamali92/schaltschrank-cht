"""
em_kopplung.py - Phase 4: EM-thermische Kopplung Python <-> Fluent.

Schleife:  T (Cu-Zonen, volumengemittelt aus CFD) -> R(T), Skin-/Proximity-Faktoren (em_modell.py)
           -> Verluste je Zone -> volumetrische Quellen in Fluent -> CFD quasi-stationaer -> neue T
           bis  max |dP|/P < EM_TOL_P  und  max |dT| < EM_TOL_T  (config.py).
Geometrie/Netz: schaltschrank_em.step / mesh/schaltschrank_em_mittel.msh.h5 (build_geometry_em.py).
Uebrige Bauteile (Umrichter, Netzteil, Schuetz) wie im Basisfall.

Aufruf:  python em_kopplung.py [--nproc 12]
Ergebnis: laeufe/em_kopplung/ (Case/Data, ergebnis.json, kopplung.csv, Ebenen-CSV),
          ergebnisse/em/kopplung_verlauf.png, logs/em_kopplung.txt
"""
import argparse
import csv
import json
import sys
import time
from pathlib import Path

import build_geometry_em as GEO
import config as C
import em_modell as EM
import solve_schaltschrank as S

HIER = Path(__file__).resolve().parent
FALL = "em_kopplung"
NETZ = HIER / "mesh" / "schaltschrank_em_mittel.msh.h5"

FESTE = ("Umrichter", "Netzteil", "Schuetz")
KERN = [f"trafo_{n.lower()}-trafo_{n.lower()}" for n, _ in C.TRAFO_TEILE_X if n.startswith("Kern")]
WICKLUNG = [f"trafo_{n.lower()}-trafo_{n.lower()}" for n, _ in C.TRAFO_TEILE_X if n.startswith("Wicklung")]
SCHIENE = {ph: f"schiene_{ph.lower()}-schiene_{ph.lower()}" for ph in C.SCHIENE_Z0}
GRUPPEN = {n: [S.ZONE[n]] for n in FESTE}
GRUPPEN.update({"Kern": KERN, "Wicklung": WICKLUNG}, **{ph: [z] for ph, z in SCHIENE.items()})


def zonenvolumen():
    vol = {}
    for name, (o, s) in GEO.koerper().items():
        vol[f"{name.lower()}-{name.lower()}"] = s[0] * s[1] * s[2]
    return vol


def stoffe_em(s):
    st = s.settings
    for name, m in (("kern", C.KERN_MAT), ("wicklung", C.WICKLUNG_MAT), ("kupfer", C.SCHIENE_MAT)):
        st.setup.materials.solid.create(name)
        mat = st.setup.materials.solid[name]
        mat.density.value = m["rho"]
        mat.specific_heat.value = m["cp"]
        mat.thermal_conductivity.value = m["k"]
    czc = st.setup.cell_zone_conditions.solid
    for z in KERN:
        czc[z].general.material = "kern"
    for z in WICKLUNG:
        czc[z].general.material = "wicklung"
    for z in SCHIENE.values():
        czc[z].general.material = "kupfer"


def quellen_setzen(s, P_zone, vol):
    czc = s.settings.setup.cell_zone_conditions.solid
    for z, P in P_zone.items():
        czc[z].sources.enable = True
        czc[z].sources.terms = {"energy": [{"option": "value", "value": P / vol[z]}]}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nproc", type=int, default=C.NPROC_SOLVER)
    a = ap.parse_args()
    ordner = HIER / "laeufe" / FALL
    ordner.mkdir(parents=True, exist_ok=True)
    (HIER / "ergebnisse" / "em").mkdir(parents=True, exist_ok=True)
    log = S.Protokoll(HIER / "logs" / f"{FALL}.txt")
    log(f"EM-thermische Kopplung ({time.strftime('%Y-%m-%d %H:%M')}), Netz {NETZ.name}")
    p = S.standard_parameter("mittel")
    p["netz"] = "em_mittel"
    p["fpsc"] = 4
    vol = zonenvolumen()
    solid_soll = {S.ZONE[n] for n in FESTE} | set(KERN) | set(WICKLUNG) | set(SCHIENE.values())
    t_start = time.time()
    s = S.starte(a.nproc, ordner)
    ergebnis = {"fall": FALL, "parameter": p, "ok": False}
    verlauf = []
    try:
        aussen = S.netz_lesen(s, NETZ, log, solid_soll=solid_soll)
        S.modelle_und_stoffe(s, p, log)
        stoffe_em(s)
        log("Feste Bauteile:")
        S.bauteile_setzen(s, p, log, nur=FESTE)
        S.numerik(s, p)
        st = s.settings
        st.solution.initialization.initialization_type = "standard"
        st.solution.initialization.defaults["temperature"] = S.K(p["T_umg"] + 10.0)
        st.solution.initialization.initialize()
        seiten = S.seiten_zuordnen(s, aussen)
        h = S.h_start(p)
        st.setup.models.radiation.model = "s2s"
        S.aussenwand_setzen(s, p, seiten, h)
        S.innenwaende_setzen(s, p, aussen, log)
        S.s2s_sichtfaktoren(s, p, ordner, log)
        r = S.reports_anlegen(s, seiten, aussen, GRUPPEN)
        st.solution.initialization.initialize()

        T = {ph: 80.0 for ph in SCHIENE} | {"wicklung": 80.0}      # Startschaetzung
        P_fest = sum(p[f"P_{n}"] for n in FESTE)
        konv_kopplung = False
        for k in range(1, C.EM_ITER_MAX + 1):
            P = EM.zonen_verluste(T)
            quellen_setzen(s, P, vol)
            Q = P_fest + sum(P.values())
            Pg = {ph: P[SCHIENE[ph]] for ph in SCHIENE} | {"Wicklung": sum(P[z] for z in WICKLUNG)}
            log(f"--- Kopplung {k}: T_ein " + ", ".join(f"{x} {v:.2f} C" for x, v in T.items())
                + " | P " + ", ".join(f"{x} {v:.2f} W" for x, v in Pg.items()) + f" | Q_ges {Q:.1f} W")
            konv, it, v, ampl, h = S.loese(s, p, r, seiten, h, log, ordner, Q=Q, monitor=f"monitor_k{k}.csv")
            T_neu = {ph: v[r["tavg"][ph]] - 273.15 for ph in SCHIENE} | {"wicklung": v[r["tavg"]["Wicklung"]] - 273.15}
            P_neu = EM.zonen_verluste(T_neu)
            dP = max(abs(P_neu[z] - P[z]) / P[z] for z in P)
            dT = max(abs(T_neu[x] - T[x]) for x in T)
            zeile = dict(kopplung=k, iterationen=it, stationaer=konv, Q_ges_W=Q,
                         bilanz_proz=100 * (Q + v["q-aussen"]) / Q, dP_rel=dP, dT_K=dT)
            zeile |= {f"P_{x}_W": val for x, val in Pg.items()}
            zeile |= {f"T_{x}_ein_C": val for x, val in T.items()}
            zeile |= {f"T_{x}_aus_C": val for x, val in T_neu.items()}
            zeile |= {f"Tmax_{g}_C": v[r["tmax"][g]] - 273.15 for g in GRUPPEN}
            zeile["Tluft_C"] = v["tavg-luft"] - 273.15
            verlauf.append(zeile)
            log(f"    -> T_aus " + ", ".join(f"{x} {val:.2f} C" for x, val in T_neu.items())
                + f" | max dP {100 * dP:.3f} %, max dT {dT:.3f} K, Bilanz {zeile['bilanz_proz']:+.3f} %")
            T = T_neu
            if dP < C.EM_TOL_P and dT < C.EM_TOL_T:
                konv_kopplung = True
                break

        with open(ordner / "kopplung.csv", "w", newline="", encoding="utf-8") as fh:
            wr = csv.DictWriter(fh, fieldnames=list(verlauf[0]))
            wr.writeheader()
            wr.writerows(verlauf)
        letzte = verlauf[-1]
        ergebnis.update(
            ok=True, kopplung_konvergiert=konv_kopplung, kopplungsschritte=len(verlauf),
            konvergiert=letzte["stationaer"], iterationen=sum(z["iterationen"] for z in verlauf),
            bilanzfehler_proz=letzte["bilanz_proz"], Q_quellen_W=letzte["Q_ges_W"],
            Tmax_C={g: letzte[f"Tmax_{g}_C"] for g in GRUPPEN}, Tluft_mittel_C=letzte["Tluft_C"],
            T_cu_C=T, P_W={x: letzte[f"P_{x}_W"] for x in list(SCHIENE) + ["Wicklung"]},
            P_kern_W=C.P_KERN, h_aussen=h)
        S.ebenen_exportieren(s, ordner, log)
        st.file.write_case_data(file_name=(ordner / f"{FALL}.cas.h5").as_posix())
        ergebnis["laufzeit_s"] = time.time() - t_start
        log("")
        log("=== ERGEBNIS EM-KOPPLUNG ===")
        log(f"Kopplung {'konvergiert' if konv_kopplung else 'NICHT konvergiert'} nach {len(verlauf)} Schritten, "
            f"{ergebnis['iterationen']} CFD-Iterationen, {ergebnis['laufzeit_s'] / 60:.1f} min")
        log(f"Energiebilanz letzter Schritt: {letzte['bilanz_proz']:+.3f} % (Q_ges {letzte['Q_ges_W']:.1f} W)")
        for g in GRUPPEN:
            log(f"  {g:10s}: Tmax {letzte[f'Tmax_{g}_C']:6.2f} C")
        log("  Cu-Temperaturen (Mittel): " + ", ".join(f"{x} {v:.2f} C" for x, v in T.items()))
        log("  Verluste: " + ", ".join(f"{x} {v:.2f} W" for x, v in ergebnis["P_W"].items())
            + f", Kern {C.P_KERN:.1f} W")
        bild(verlauf)
    except Exception as e:
        import traceback
        log("!!! FEHLER: " + repr(e))
        log(traceback.format_exc())
        ergebnis["fehler"] = repr(e)
    finally:
        (ordner / "ergebnis.json").write_text(json.dumps(ergebnis, indent=2, default=float), encoding="utf-8")
        try:
            s.exit()
        except Exception:
            pass


def bild(verlauf):
    import plotstil as PS
    PS.stil()
    plt = PS.plt
    k = [int(z["kopplung"]) for z in verlauf]
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(8.8, 3.9))
    farben = {"L1": PS.KATEGORIE[0], "L2": PS.KATEGORIE[1], "L3": PS.KATEGORIE[2], "Wicklung": PS.KATEGORIE[3]}
    # Verluste nach dem letzten Schritt (mit T_aus) ergaenzen -> Verlauf 1..n+1
    for x, fa in farben.items():
        tx = "wicklung" if x == "Wicklung" else x
        P = [float(z[f"P_{x}_W"]) for z in verlauf]
        a1.plot(k, [100 * (v / P[0] - 1) for v in P], "o-", color=fa, ms=5, label=x)
        T = [float(verlauf[0][f"T_{tx}_ein_C"])] + [float(z[f"T_{tx}_aus_C"]) for z in verlauf]
        a2.plot([0] + k, T, "o-", color=fa, ms=5, label=x)
    a1.axhline(0, color=PS.TEXT2, lw=0.8)
    a1.set_xlabel("Kopplungsschritt")
    a1.set_ylabel("Verlustaenderung gegenueber Schritt 1 in %")
    a1.set_title("Verluste R(T), Skin/Proximity (Eingang je Schritt)", fontsize=9)
    a2.set_xlabel("Kopplungsschritt (0 = Startwert)")
    a2.set_ylabel("mittlere Cu-Temperatur in °C")
    a2.set_title("CFD-Temperatur (Ergebnis je Schritt)", fontsize=9)
    a1.set_xticks(k)
    a2.set_xticks([0] + k)
    a2.legend(fontsize=8)
    fig.suptitle("EM-thermische Kopplung Python <-> Fluent (Startwert 80 °C); L1 und L3 liegen symmetrisch und decken sich", fontsize=9.5)
    fig.savefig(HIER / "ergebnisse" / "em" / "kopplung_verlauf.png")
    plt.close(fig)


if __name__ == "__main__":
    if "--nur-bild" in sys.argv:
        with open(HIER / "laeufe" / FALL / "kopplung.csv", encoding="utf-8") as fh:
            bild(list(csv.DictReader(fh)))
    else:
        main()
