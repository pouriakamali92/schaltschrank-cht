"""
check_env.py - Umgebungscheck fuer das Schaltschrank-CFD-Projekt (Windows)

Reines Lesen/Pruefen: aendert nichts am System, startet keine Simulation.
Aufruf:  python check_env.py
Ausgabe wird zusaetzlich in check_env_report.txt gespeichert.
"""

import os
import platform
import shutil
import subprocess
import sys
from datetime import datetime

lines = []


def out(text=""):
    print(text)
    lines.append(text)


def run(cmd, timeout=20):
    try:
        r = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, shell=False)
        return (r.stdout or r.stderr).strip()
    except Exception as e:
        return f"<Fehler: {e}>"


def section(title):
    out()
    out("=" * 60)
    out(title)
    out("=" * 60)


# ---------------------------------------------------------------- System
section("1. System")
out(f"Datum/Zeit     : {datetime.now():%Y-%m-%d %H:%M:%S}")
out(f"Betriebssystem : {platform.platform()}")
out(f"Rechnername    : {platform.node()}")
out(f"CPU-Kerne      : {os.cpu_count()}")
out(f"Arbeitsordner  : {os.getcwd()}")

try:
    import ctypes

    class MEMSTAT(ctypes.Structure):
        _fields_ = [
            ("dwLength", ctypes.c_ulong), ("dwMemoryLoad", ctypes.c_ulong),
            ("ullTotalPhys", ctypes.c_ulonglong), ("ullAvailPhys", ctypes.c_ulonglong),
            ("ullTotalPageFile", ctypes.c_ulonglong), ("ullAvailPageFile", ctypes.c_ulonglong),
            ("ullTotalVirtual", ctypes.c_ulonglong), ("ullAvailVirtual", ctypes.c_ulonglong),
            ("sullAvailExtendedVirtual", ctypes.c_ulonglong),
        ]

    m = MEMSTAT()
    m.dwLength = ctypes.sizeof(MEMSTAT)
    ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(m))
    out(f"RAM gesamt     : {m.ullTotalPhys / 1024**3:.1f} GB")
    out(f"RAM frei       : {m.ullAvailPhys / 1024**3:.1f} GB")
except Exception:
    out("RAM            : nicht ermittelbar")

total, used, free = shutil.disk_usage(os.getcwd())
out(f"Platte (frei)  : {free / 1024**3:.1f} GB von {total / 1024**3:.1f} GB")

# ---------------------------------------------------------------- Python
section("2. Python")
out(f"Python-Version : {sys.version.split()[0]}")
out(f"Interpreter    : {sys.executable}")

# ---------------------------------------------------------------- Pakete
section("3. Python-Pakete")
packages = [
    ("ansys.fluent.core", "ansys-fluent-core (PyFluent)"),
    ("ansys.fluent.visualization", "ansys-fluent-visualization"),
    ("numpy", "numpy"),
    ("scipy", "scipy"),
    ("pandas", "pandas"),
    ("matplotlib", "matplotlib"),
]
import importlib

for mod, label in packages:
    try:
        m_ = importlib.import_module(mod)
        out(f"[OK]    {label:38s} {getattr(m_, '__version__', '')}")
    except Exception:
        out(f"[FEHLT] {label}")

# ---------------------------------------------------------------- Fluent
section("4. Ansys Fluent")
awp = {k: v for k, v in os.environ.items() if k.upper().startswith("AWP_ROOT")}
if awp:
    for k, v in sorted(awp.items()):
        out(f"{k} = {v}")
        fl = os.path.join(v, "fluent", "ntbin", "win64", "fluent.exe")
        out(f"   fluent.exe vorhanden: {os.path.exists(fl)}")
else:
    out("Keine AWP_ROOT-Umgebungsvariable gefunden.")

fluent_path = shutil.which("fluent")
out(f"'fluent' im PATH: {fluent_path if fluent_path else 'nein'}")

lic = {k: v for k, v in os.environ.items() if "LICENSE" in k.upper() or "LMD" in k.upper()}
for k, v in lic.items():
    out(f"{k} = {v}")

# ---------------------------------------------------------------- Werkzeuge
section("5. Weitere Werkzeuge")
for tool in ["git", "node", "npm", "claude", "pip"]:
    p = shutil.which(tool)
    if p:
        ver = run([p, "--version"])
        out(f"[OK]    {tool:7s} {ver.splitlines()[0] if ver else ''}")
    else:
        out(f"[FEHLT] {tool}")

# ---------------------------------------------------------------- Schluss
section("6. Hinweis")
out("Bitte den Inhalt von check_env_report.txt komplett in den Chat kopieren.")
out("Keine Lizenzserver-Adressen oder interne Pfade posten, falls diese")
out("vertraulich sind - einfach durch XXX ersetzen.")

with open("check_env_report.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(lines))
