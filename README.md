# threeBodySim
## Python 3.12

## Drei-Körper-Simulation (by Timo Streich)

Diese Python-Anwendung simuliert das Dreikörperproblem in drei Dimensionen und stellt die Bewegung der Körper als Animation dar. Die Körper interagieren über die Gravitation, und die Bahnen werden numerisch berechnet.
Voraussetzungen

Stellen Sie sicher, dass die folgenden Python-Bibliotheken installiert sind:

    numpy

    scipy

    matplotlib

Sie können diese mit dem folgenden Befehl installieren:

    pip install numpy scipy matplotlib

### Verwendung

    Simulation ausführen:

    Starten Sie das Skript index.py mit Python:
    python index.py

    Folgende Argumente sind vorhanden:
    -help/-h: zeigt die Kommands an
    -line: eine linie zu den Punkten wird angezeigt
    -shortline: eine linie mit der Länge MAX_LENGTH wird zu den Punkten angezeigt
    -slow: die Simulation läuft langsamer
    -planet: ein Planet wird hinzugefügt


    Es öffnet sich ein 3D Fenster, inwelchem man sich mit der Maus bewegen kann.

### Funktionsweise

    Initialisierung:

    Die Anfangsbedingungen der drei Körper (Positionen und Geschwindigkeiten) sind im Code definiert.

    Numerische Integration:

    Die Bewegungsgleichungen werden mit der Funktion odeint aus scipy.integrate löst, um die Positionen der Körper über die Zeit zu berechnen.

    Visualisierung:

    Die berechneten Positionen werden mit matplotlib visualisiert. Die Animation zeigt sowohl die aktuellen Positionen der Körper als auch die Bahnen, die sie im Laufe der Zeit beschreiben.


Sie können die Simulation anpassen, indem Sie die Anfangsbedingungen oder die Massen der Körper ändern. Suchen Sie nach der Funktion initial_conditions im Skript, um Positionen und Geschwindigkeiten zu modifizieren:
### Beispiel: Positionen
    r1 = [-1.0, 0.0, 0.0]
    r2 = [1.0, 0.0, 0.0]
    r3 = [0.0, 1.0, 0.0]

### Beispiel: Geschwindigkeiten
    v1 = [0.0, -0.5, 0.0]
    v2 = [0.0, 0.5, 0.0]
    v3 = [0.0, 0.0, 0.5]

## Lizenz

Dieses Projekt steht unter der MIT-Lizenz.



