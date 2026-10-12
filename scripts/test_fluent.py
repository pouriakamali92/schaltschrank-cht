"""
test_fluent.py - Minimaltest: startet Fluent headless ueber PyFluent und beendet es wieder.
Aufruf:  python test_fluent.py
"""

import ansys.fluent.core as pyfluent

print("PyFluent:", pyfluent.__version__)

solver = pyfluent.launch_fluent(
    mode="solver",
    dimension=3,
    precision="double",
    processor_count=4,   # bei Lizenzfehler auf 2 oder 1 reduzieren
    ui_mode="no_gui",
    cleanup_on_exit=True,
)

print("Fluent laeuft. Version:", solver.get_fluent_version())

solver.exit()
print("Test erfolgreich.")
