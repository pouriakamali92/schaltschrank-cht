"""
run_mesh.py - fuehrt die komplette Netzphase automatisch aus:
  grob -> Zonenpruefung -> mittel -> Zonenpruefung -> fein -> Zonenpruefung
Bei einem Fehler wird abgebrochen. Am Ende steht alles Wichtige in
    logs\\zusammenfassung.txt     (nur diese Datei in den Chat kopieren)

Aufruf:  D:\\Pouria\\python\\python.exe run_mesh.py [grob mittel fein]
"""
import os
import re
import subprocess
import sys
import time
from pathlib import Path

HIER = Path(__file__).resolve().parent
LOGS = HIER / "logs"
LOGS.mkdir(exist_ok=True)
stufen = sys.argv[1:] or ["grob", "mittel", "fein"]
env = dict(os.environ, PYTHONIOENCODING="utf-8")

zusammen = []


def z(text=""):
    print(text, flush=True)
    zusammen.append(text)


def lauf(skript, stufe, logname):
    log = LOGS / logname
    t0 = time.time()
    with open(log, "w", encoding="utf-8", errors="replace") as f:
        try:
            r = subprocess.run([sys.executable, "-u", str(HIER / skript), stufe],
                               cwd=HIER, stdout=f, stderr=subprocess.STDOUT, env=env, timeout=3600)
            code = r.returncode
        except subprocess.TimeoutExpired:
            f.write("\n!!! TIMEOUT nach 3600 s\n")
            code = 99
    return code, log.read_text(encoding="utf-8", errors="replace"), time.time() - t0


def zahl(muster, text, letzte=True):
    m = re.findall(muster, text)
    return (m[-1] if letzte else m[0]) if m else "?"


def auszug(text, n=45):
    return "\n".join(text.strip().splitlines()[-n:])


z(f"Netzphase gestartet, Stufen: {stufen}")
fehler = False
for st in stufen:
    z(f"\n===== {st} =====")
    code, text, dauer = lauf("mesh_schaltschrank.py", st, f"mesh_{st}.log")
    if code != 0:
        z(f"FEHLER beim Vernetzen ({dauer:.0f} s). Letzte Zeilen des Logs:")
        z(auszug(text))
        fehler = True
        break
    z(f"Vernetzen ok ({dauer:.0f} s)")
    for zeile in re.findall(r"^\s*Regionen \((?:Vorgabe|gesetzt|danach)\):.*$", text, re.M):
        z("  " + zeile.strip())
    z("  Zellen             : " + zahl(r"(\d+) cells were created", text))
    skew = zahl(r"quality limits \(min max ave\) = \(\S+ ([0-9.]+)", text)
    oq = zahl(r"Minimum Orthogonal Quality =\s*([0-9.eE+-]+)", text)
    z("  Max. Skewness      : " + skew + "  (Endkontrolle alle Zellen)")
    z("  Min. Orthog.       : " + oq)
    try:
        gate = float(skew) <= 0.85 and float(oq) > 0.15
    except ValueError:
        gate = False
    z("  GATE Qualitaet     : " + ("OK" if gate else "VERLETZT (Skewness <= 0.85, OQ > 0.15)"))
    z("  Max. Aspect Ratio  : " + zahl(r"Maximum Aspect Ratio =\s*([0-9.eE+-]+)", text))
    tab = re.findall(r"^\s*\d+\s+\S+-\S+\s+(?:fluid|solid)\s+\d+.*$", text, re.M)
    if tab:
        z("  Zellzonen (id name typ anzahl tet hex pyr wedge poly poly-prism hex-prism):")
        for zeile in tab[-5:]:
            z("    " + " ".join(zeile.split()))
    pr = re.findall(r"PRISMENSCHICHTEN IM NETZ:.*", text)
    z("  " + (pr[-1].strip() if pr else "PRISMENSCHICHTEN: keine Angabe im Log"))
    for h in re.findall(r"^\s*>>> .*Prismen fehlgeschlagen.*$", text, re.M):
        z("  " + h.strip())

    code, text2, dauer = lauf("pruefe_netz.py", st, f"pruefung_{st}.log")
    zeilen = [l for l in text2.splitlines() if l.startswith(("ZONEN", "ZONENPRUEFUNG"))]
    for l in zeilen:
        z("  " + l)
    if code != 0:
        z("FEHLER bei der Zonenpruefung. Letzte Zeilen des Logs:")
        z(auszug(text2, 25))
        fehler = True
        break

z("\n" + ("ABGEBROCHEN (siehe Fehler oben)." if fehler else "ALLES OK: alle Netze erzeugt und geprueft."))
(LOGS / "zusammenfassung.txt").write_text("\n".join(zusammen), encoding="utf-8")
print(f"\nGeschrieben: {LOGS / 'zusammenfassung.txt'}")
sys.exit(1 if fehler else 0)
