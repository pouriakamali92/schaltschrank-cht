# Wärme im Schaltschrank: CHT-, Strahlungs- und elektrisch-thermische Simulation

Strömungs-, Strahlungs- und elektrisch-thermische Simulation eines Schaltschranks mit Ansys Fluent, mit Gegenprüfungen, Szenarien und einem schnellen Ersatzmodell.

**Art:** persönliches Selbststudienprojekt, 2026
**Software:** Ansys Fluent 2024 R2 (PyFluent), Python, MATLAB/Simulink

> **Modellstudie ohne Messdaten. Verifiziert, nicht validiert.** Verluste, Wärmeleitfähigkeiten und Oberflächeneigenschaften sind dokumentierte, plausible Annahmen. Die Ergebnisse zeigen, wie ein solcher Schrank grundsätzlich reagiert, nicht wie ein bestimmter realer Schrank.

## Aufgabe

Schrank 0,8 × 1,2 × 0,4 m mit Umrichter, Netzteil, Schütz und Trafo, zusammen 250 W Verlustleistung, ohne Lüfter. Leitfragen: Wie heiß werden die Bauteile bei 35 °C Umgebung? Wie sicher ist die Rechnung? Was bringen Lüfter mit Filter, anderes Layout und 45 °C Umgebung? Wie ändern sich elektrische Verluste mit Frequenz und Temperatur?

## Vorgehen

- Geometrie, Verluste und Materialwerte in einer Konfigurationsdatei, Annahmen dokumentiert
- Gitter mit Prismenschichten in drei Feinheiten (45.000, 135.000, 353.000 Zellen), automatisch erzeugt und qualitätsgeprüft
- Gekoppelte Simulation (CHT) in Fluent: Konvektion, Wärmeleitung in den Bauteilen, Strahlung zwischen allen Flächen, vollständig per Skript gesteuert
- Gegenprüfungen: Energiebilanz je Rechnung, Gitterstudie, Annahmen-Empfindlichkeit, Modellvarianten, Handrechnung und unabhängiges Handmodell, Plausibilitäts-Audit mit 77 Einzelprüfungen (keine Fehler, 9 erklärte Hinweise)
- Elektrisch-thermische Kopplung für Stromschienen und Trafo-Wicklung (Skin- und Proximity-Effekt, analytisch), iterativ mit der Wärmerechnung abgeglichen
- Ersatzmodell als kleines Wärmenetz aus den Fluent-Rechnungen, in Python und Simulink
- Szenarien: Lüfter mit Kennlinie und Filter, Layout B, 45 °C Umgebung

## Ergebnisse (Auswahl)

**Basisfall:** mittlere Lufttemperatur 50 °C, Außenwand 43,5 °C, etwa 13 K Unterschied oben/unten. Etwa zwei Drittel der Wärme verlassen den Schrank außen als Strahlung.

| Bauteil (heißester Punkt) | CFD | Handmodell (Bandbreite) |
|---|---|---|
| Umrichter | 78,0 °C | 71,6 – 83,5 °C |
| Netzteil | 87,2 °C | 78,1 – 93,5 °C |
| Schütz | 103,8 °C | 91,3 – 114,6 °C |
| Trafo | 76,3 °C | 73,3 – 86,3 °C |

Die Bänder sind 12 bis 23 K breit: Der Test findet grobe Fehler, keine 2 bis 3 K.

| Szenario | Schütz | Luft |
|---|---|---|
| Basis | 103,8 °C | 50,0 °C |
| 45 °C Umgebung | 111,7 °C | 59,3 °C |
| Layout B | 100,8 °C | 49,7 °C |
| Filterlüfter | 86,0 °C | 39,8 °C |
| Filterlüfter + 45 °C | 94,8 °C | 49,8 °C |

- Filterlüfter (rund 67 m³/h) senkt die Lufttemperatur um 10 K.
- Der Emissionsgrad ist der stärkste Hebel, weil Strahlung dominiert.
- Elektrisch-thermisch: Eindringtiefe in Kupfer bei 50 Hz etwa 9,3 mm, bei 1 kHz 2,1 mm. Oberschwingungen einer 6-pulsigen Last erhöhen die Schienenverluste um etwa 8 %.
- Ersatzmodell: 1,14 K Abweichung im Prüffall, der nicht zur Anpassung diente. Es braucht etwa 0,1 s je Fall (Simulink) gegenüber 5 bis 11 min stationär in Fluent. Der Vorteil gilt pro Auswertung nach dem Aufbau, denn das Modell wurde aus den Fluent-Rechnungen angepasst.

## Belastbarkeit und Grenzen

- Energiebilanz höchstens 0,43 % in allen 23 Rechnungen. Ein anfänglicher Fehler von 1,9 % (Standardeinstellung des Strahlungslösers) wurde gefunden und behoben, danach wurden alle Läufe neu gerechnet.
- Gitter: Temperaturen ändern sich zwischen den drei Netzen um höchstens 0,75 K. Ein formaler GCI ist nur für den Trafo sauber bestimmbar (etwa 0,16 K).
- Keine Messdaten. Annahmen: Verluste 120 / 40 / 25 / 65 W, Bauteil-Wärmeleitfähigkeit 20 / 5 / 1,5 / 10 W/(m·K), Emissionsgrad 0,9.
- Bauteile sind homogene Blöcke. Der heißeste Punkt ist ein Blockwert, keine Chip- oder Kontakttemperatur.
- Strömung im Übergang laminar/turbulent (Ra ≈ 10⁹), die Lösung schwankt leicht.
- Elektrik nur analytisch, das 1D-Schienenmodell ist nur bis etwa 500 Hz belastbar.

Portfolio mit Abbildungen: <https://pouria-kamali.de/projekt-schaltschrank.php>

## Abbildungen

- [basis_monitore.png](images/basis_monitore.png)
- [basis_schichtung.png](images/basis_schichtung.png)
- [basis_temperatur_x200.png](images/basis_temperatur_x200.png)
- [basis_temperatur_x650.png](images/basis_temperatur_x650.png)
- [basis_temperatur_z060.png](images/basis_temperatur_z060.png)
- [basis_temperatur_z300.png](images/basis_temperatur_z300.png)
- [em_em_kopplung_temperatur_x650.png](images/em_em_kopplung_temperatur_x650.png)
- [em_em_kopplung_temperatur_z060.png](images/em_em_kopplung_temperatur_z060.png)
- [em_kopplung_verlauf.png](images/em_kopplung_verlauf.png)
- [em_schienen_rac_rdc.png](images/em_schienen_rac_rdc.png)
- [em_trafo_dowell.png](images/em_trafo_dowell.png)
- [rom_rom_fahne.png](images/rom_rom_fahne.png)
- [rom_rom_simulink.png](images/rom_rom_simulink.png)
- [rom_rom_stationaer.png](images/rom_rom_stationaer.png)
- [rom_rom_transient.png](images/rom_rom_transient.png)
- [szenarien_szenarien_tmax.png](images/szenarien_szenarien_tmax.png)
- [szenarien_sz_layoutB_temperatur_z060.png](images/szenarien_sz_layoutB_temperatur_z060.png)
- [szenarien_sz_luefter_temperatur_z060.png](images/szenarien_sz_luefter_temperatur_z060.png)
- [szenarien_sz_luefter_temperatur_z300.png](images/szenarien_sz_luefter_temperatur_z300.png)
- [vv_gci.png](images/vv_gci.png)
- [vv_vv_sensitivitaet.png](images/vv_vv_sensitivitaet.png)

## Skripte

- [audit_fluent_lesen.py](scripts/audit_fluent_lesen.py)
- [build_geometry.py](scripts/build_geometry.py)
- [build_geometry_em.py](scripts/build_geometry_em.py)
- [check_env.py](scripts/check_env.py)
- [config.py](scripts/config.py)
- [config_alt.py](scripts/config_alt.py)
- [em_kopplung.py](scripts/em_kopplung.py)
- [em_modell.py](scripts/em_modell.py)
- [gci.py](scripts/gci.py)
- [handrechnung.py](scripts/handrechnung.py)
- [mesh_schaltschrank.py](scripts/mesh_schaltschrank.py)
- [plausibilitaet_audit.py](scripts/plausibilitaet_audit.py)
- [plotstil.py](scripts/plotstil.py)
- [plots_fall.py](scripts/plots_fall.py)
- [pruefe_netz.py](scripts/pruefe_netz.py)
- [rom.py](scripts/rom.py)
- [rom_export_erzeugen.py](scripts/rom_export_erzeugen.py)
- [rom_fahne.py](scripts/rom_fahne.py)
- [rom_laeufe.py](scripts/rom_laeufe.py)
- [rom_matlab.m](scripts/rom_matlab.m)
- [rom_simulink.m](scripts/rom_simulink.m)
- [rom_transient.py](scripts/rom_transient.py)
- [rom_vergleich.py](scripts/rom_vergleich.py)
- [run_mesh.py](scripts/run_mesh.py)
- [sammle_ergebnisse.py](scripts/sammle_ergebnisse.py)
- [solver_info.py](scripts/solver_info.py)
- [solve_schaltschrank.py](scripts/solve_schaltschrank.py)
- [szenarien.py](scripts/szenarien.py)
- [test_fluent.py](scripts/test_fluent.py)
- [vv_laeufe.py](scripts/vv_laeufe.py)
