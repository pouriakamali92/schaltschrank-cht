"""
rom_transient.py - Phase 5: transienter Fluent-Lauf fuer die ROM-Pruefung (Netz grob).

Anfangszustand: alles auf Umgebungstemperatur, Luft in Ruhe. Bei t = 0 werden alle Verlustleistungen
(Basisfall) eingeschaltet. Aussen-h fest auf den nachgefuehrten Werten des stationaeren Basislaufs
(basis_grob) - Vereinfachung, dokumentiert.
Zeitschritt 30 s (implizit, 1. Ordnung), max. 12 Iterationen je Schritt, Dauer 8 h.

Aufruf: python rom_transient.py [--stunden 8] [--dt 30]
Ergebnis: laeufe/rom_transient/transient.csv, logs/rom_transient.txt
"""
import argparse
import json
import time
from pathlib import Path

import config as C
import solve_schaltschrank as S

HIER = Path(__file__).resolve().parent
FALL = "rom_transient"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--stunden", type=float, default=8.0)
    ap.add_argument("--dt", type=float, default=30.0)
    ap.add_argument("--nproc", type=int, default=C.NPROC_SOLVER)
    a = ap.parse_args()
    ordner = HIER / "laeufe" / FALL
    ordner.mkdir(parents=True, exist_ok=True)
    log = S.Protokoll(HIER / "logs" / f"{FALL}.txt")
    p = S.standard_parameter("grob")
    basis = json.loads((HIER / "laeufe" / "basis_grob" / "ergebnis.json").read_text(encoding="utf-8"))
    h = basis["h_aussen"]
    log(f"Transient ({time.strftime('%Y-%m-%d %H:%M')}): Netz grob, dt {a.dt} s, {a.stunden} h, h aussen fest {h}")
    t0 = time.time()
    s = S.starte(a.nproc, ordner)
    try:
        aussen = S.netz_lesen(s, HIER / "mesh" / "schaltschrank_grob.msh.h5", log)
        S.modelle_und_stoffe(s, p, log)
        S.bauteile_setzen(s, p, log)
        S.numerik(s, p)
        st = s.settings
        st.solution.initialization.initialization_type = "standard"
        st.solution.initialization.defaults["temperature"] = S.K(p["T_umg"])
        st.solution.initialization.initialize()
        seiten = S.seiten_zuordnen(s, aussen)
        st.setup.models.radiation.model = "s2s"
        S.aussenwand_setzen(s, p, seiten, h)
        S.innenwaende_setzen(s, p, aussen, log)
        S.s2s_sichtfaktoren(s, p, ordner, log)
        r = S.reports_anlegen(s, seiten, aussen)
        namen = S.alle_namen(r)
        # transient
        st.setup.general.solver.time = "unsteady-1st-order"
        tc = st.solution.run_calculation.transient_controls
        tc.time_step_size = a.dt
        st.solution.initialization.initialize()
        n_ges = int(round(a.stunden * 3600 / a.dt))
        block = 10                          # Ausgabe alle 10 Zeitschritte
        zeilen = []
        v = S.berechne(st, namen)
        zeilen.append(dict(v, t=0.0))
        n = 0
        while n < n_ges:
            k = min(block, n_ges - n)
            st.solution.run_calculation.dual_time_iterate(time_step_count=k, max_iter_per_step=12)
            n += k
            v = S.berechne(st, namen)
            zeilen.append(dict(v, t=n * a.dt))
            if n % 120 == 0 or n == n_ges:
                log(f"  t = {n * a.dt / 3600:5.2f} h: " + " ".join(
                    f"{b[:4]} {v[r['tavg'][b]] - 273.15:6.2f}" for b in S.BAUTEILE)
                    + f" | Luft {v['tavg-luft'] - 273.15:6.2f} | Q_aussen {-v['q-aussen']:7.1f} W"
                    + f" | {time.time() - t0:6.0f} s")
        spalten = ["t"] + namen
        with open(ordner / "transient.csv", "w", encoding="utf-8") as fh:
            fh.write(",".join(spalten) + "\n")
            for z in zeilen:
                fh.write(",".join(f"{z.get(c, float('nan')):.6g}" for c in spalten) + "\n")
        st.file.write_case_data(file_name=(ordner / f"{FALL}.cas.h5").as_posix())
        laufzeit = time.time() - t0
        (ordner / "ergebnis.json").write_text(json.dumps(dict(ok=True, laufzeit_s=laufzeit, dt=a.dt,
                                                              stunden=a.stunden, h_aussen=h), indent=2))
        log(f"Fertig: Laufzeit {laufzeit / 60:.1f} min fuer {a.stunden} h physikalische Zeit")
    except Exception as e:
        import traceback
        log("!!! FEHLER: " + repr(e))
        log(traceback.format_exc())
    finally:
        try:
            s.exit()
        except Exception:
            pass


if __name__ == "__main__":
    main()
