import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
import argparse


MAX_LENGTH = 750  # Hier die maximale Länge der Linie festlegen

# Gravitationskonstante
G = 1.0

# Massen der drei Körper und des Planeten
m1, m2, m3 = 1.0, 1.0, 1.0
m_planet = 0.01  # Geringe Masse des Planeten

# Anfangsbedingungen: Positionen und Geschwindigkeiten der drei Körper und des Planeten
def initial_conditions(planet):
    # Positionen
    r1 = [-5.0, 0.0, 0.0]
    r2 = [5.0, -5.0, 0.0]
    r3 = [0.0, 5.0, 0.0]

    # Geschwindigkeiten
    v1 = [-0.1, 0.0, 0.4]
    v2 = [0.0, 0.1, -0.2]
    v3 = [0.1, -0.1, -0.2]

    if planet:
        # Der Planet umkreist Körper 1
        r_planet = [5.2, 0.0, 0.0]  # Position leicht außerhalb von Körper 1
        v_planet = [0.0, 0.5, 0.0]  # Geschwindigkeit für eine kreisförmige Bahn
        return np.array(r1 + r2 + r3 + r_planet + v1 + v2 + v3 + v_planet)
    else:
        return np.array(r1 + r2 + r3 + v1 + v2 + v3)

# Funktion, die die Bewegungsgleichungen beschreibt
def derivatives(state, t, G, m1, m2, m3, m_planet, planet):
    # Entpacke den Zustand
    if planet:
        x1, y1, z1, x2, y2, z2, x3, y3, z3, x_planet, y_planet, z_planet, \
        vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3, vx_planet, vy_planet, vz_planet = state

        # Positionen als Arrays
        r1 = np.array([x1, y1, z1])
        r2 = np.array([x2, y2, z2])
        r3 = np.array([x3, y3, z3])
        r_planet = np.array([x_planet, y_planet, z_planet])

        # Abstandsvektoren
        r12 = np.linalg.norm(r2 - r1)
        r13 = np.linalg.norm(r3 - r1)
        r23 = np.linalg.norm(r3 - r2)
        r1p = np.linalg.norm(r_planet - r1)
        r2p = np.linalg.norm(r_planet - r2)
        r3p = np.linalg.norm(r_planet - r3)

        # Gravitationsbeschleunigungen
        a1 = G * m2 * (r2 - r1) / r12 ** 3 + G * m3 * (r3 - r1) / r13 ** 3 + G * m_planet * (r_planet - r1) / r1p ** 3
        a2 = G * m1 * (r1 - r2) / r12 ** 3 + G * m3 * (r3 - r2) / r23 ** 3 + G * m_planet * (r_planet - r2) / r2p ** 3
        a3 = G * m1 * (r1 - r3) / r13 ** 3 + G * m2 * (r2 - r3) / r23 ** 3 + G * m_planet * (r_planet - r3) / r3p ** 3
        a_planet = G * m1 * (r1 - r_planet) / r1p ** 3 + G * m2 * (r2 - r_planet) / r2p ** 3 + G * m3 * (r3 - r_planet) / r3p ** 3

        # Rückgabe der Ableitungen
        return np.array([
            vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3, vx_planet, vy_planet, vz_planet,
            a1[0], a1[1], a1[2], a2[0], a2[1], a2[2], a3[0], a3[1], a3[2], a_planet[0], a_planet[1], a_planet[2]
        ])
    else:
        x1, y1, z1, x2, y2, z2, x3, y3, z3, vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3 = state

        # Positionen als Arrays
        r1 = np.array([x1, y1, z1])
        r2 = np.array([x2, y2, z2])
        r3 = np.array([x3, y3, z3])

        # Abstandsvektoren
        r12 = np.linalg.norm(r2 - r1)
        r13 = np.linalg.norm(r3 - r1)
        r23 = np.linalg.norm(r3 - r2)

        # Gravitationsbeschleunigungen
        a1 = G * m2 * (r2 - r1) / r12 ** 3 + G * m3 * (r3 - r1) / r13 ** 3
        a2 = G * m1 * (r1 - r2) / r12 ** 3 + G * m3 * (r3 - r2) / r23 ** 3
        a3 = G * m1 * (r1 - r3) / r13 ** 3 + G * m2 * (r2 - r3) / r23 ** 3

        # Rückgabe der Ableitungen
        return np.array([
            vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3,
            a1[0], a1[1], a1[2], a2[0], a2[1], a2[2], a3[0], a3[1], a3[2]
        ])

# Argumentparser zur Steuerung der Simulation
parser = argparse.ArgumentParser(description="3-Körper-Simulation")
parser.add_argument("-line", action="store_true", help="Zeige die Linien der Bahnen an")
parser.add_argument("-shortline", action="store_true", help="Zeige die Linien der Bahnen mit der Länge MAX_LENGTH an")
parser.add_argument("-slow", action="store_true", help="Verlangsamt die Animation")
parser.add_argument("-planet", action="store_true", help="Füge einen Planeten hinzu")
parser.add_argument("-help", action="store_true", help="Zeigt Hilfeinformationen an")
args = parser.parse_args()

if args.help:
    parser.print_help()
    exit()

# Anfangsbedingungen laden
initial_state = initial_conditions(args.planet)

# Visualisierung der Bahnen mit Animation
size = 10
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.set_xlim([-size, size])
ax.set_ylim([-size, size])
ax.set_zlim([-size, size])
ax.set_title("Live-Simulation des Dreikörperproblems")
ax.set_xlabel("x-Position")
ax.set_ylabel("y-Position")
ax.set_zlabel("z-Position")

# Linien und Punkte für die Körper
line_planet = None  # Initialisierung, um Referenzprobleme zu vermeiden
if args.line:
    line1, = ax.plot([], [], [], lw=1, label="Sonne 1", color="blue")
    line2, = ax.plot([], [], [], lw=1, label="Sonne 2", color="green")
    line3, = ax.plot([], [], [], lw=1, label="Sonne 3", color="orange")
    if args.planet:
        line_planet, = ax.plot([], [], [], lw=1, label="Planet", color="red")

if args.shortline:
    line1, = ax.plot([], [], [], lw=1, label="Sonne 1", color="blue")
    line2, = ax.plot([], [], [], lw=1, label="Sonne 2", color="green")
    line3, = ax.plot([], [], [], lw=1, label="Sonne 3", color="orange")
    if args.planet:
        line_planet, = ax.plot([], [], [], lw=1, label="Planet", color="red")


else:
    line1 = line2 = line3 = line_planet = None

point1, = ax.plot([], [], [], 'o', color="blue")
point2, = ax.plot([], [], [], 'o', color="green")
point3, = ax.plot([], [], [], 'o', color="orange")
if args.planet:
    point_planet, = ax.plot([], [], [], 'o', color="red")
else:
    point_planet = None  # Sicherstellen, dass es definiert ist

# Speicher für die Bahnen
trail1, trail2, trail3, trail_planet = [], [], [], []
solution = [initial_state]

# Initialisierungsfunktion
def init():
    if args.line:
        line1.set_data([], [])
        line1.set_3d_properties([])
        line2.set_data([], [])
        line2.set_3d_properties([])
        line3.set_data([], [])
        line3.set_3d_properties([])
        if args.planet:
            line_planet.set_data([], [])
            line_planet.set_3d_properties([])

    if args.shortline:
        line1.set_data([], [])
        line1.set_3d_properties([])
        line2.set_data([], [])
        line2.set_3d_properties([])
        line3.set_data([], [])
        line3.set_3d_properties([])
        if args.planet:
            line_planet.set_data([], [])
            line_planet.set_3d_properties([])

    point1.set_data([], [])
    point1.set_3d_properties([])
    point2.set_data([], [])
    point2.set_3d_properties([])
    point3.set_data([], [])
    point3.set_3d_properties([])
    if args.planet:
        point_planet.set_data([], [])
        point_planet.set_3d_properties([])
    return (line1, line2, line3, line_planet, point1, point2, point3, point_planet) if args.planet and args.line else (point1, point2, point3, point_planet)

# Update-Funktion für die Animation
def update(frame):
    global trail1, trail2, trail3, trail_planet, solution

    # Berechne neue Lösungen, falls nötig
    if frame >= len(solution):
        t = np.linspace(0, 1000, 5000)  # Kurzer Zeitabschnitt
        new_solution = odeint(derivatives, solution[-1], t, args=(G, m1, m2, m3, m_planet, args.planet))
        solution.extend(new_solution[1:])  # Vermeide Dopplung des letzten Werts

    # Hole aktuelle Positionen
    state = solution[frame]
    if args.planet:
        x1, y1, z1, x2, y2, z2, x3, y3, z3, x_planet, y_planet, z_planet = state[:12]
        trail_planet.append((x_planet, y_planet, z_planet))
    else:
        x1, y1, z1, x2, y2, z2, x3, y3, z3 = state[:9]

    trail1.append((x1, y1, z1))
    trail2.append((x2, y2, z2))
    trail3.append((x3, y3, z3))

    # Update der Linien
    if args.line:
        line1.set_data([pos[0] for pos in trail1], [pos[1] for pos in trail1])
        line1.set_3d_properties([pos[2] for pos in trail1])
        line2.set_data([pos[0] for pos in trail2], [pos[1] for pos in trail2])
        line2.set_3d_properties([pos[2] for pos in trail2])
        line3.set_data([pos[0] for pos in trail3], [pos[1] for pos in trail3])
        line3.set_3d_properties([pos[2] for pos in trail3])
        if args.planet:
            line_planet.set_data([pos[0] for pos in trail_planet], [pos[1] for pos in trail_planet])
            line_planet.set_3d_properties([pos[2] for pos in trail_planet])

    if args.shortline:
        # Linie 1
        if len(trail1) > MAX_LENGTH:
            trail1 = trail1[-MAX_LENGTH:]  # Behalte nur die letzten MAX_LENGTH Punkte
        line1.set_data([pos[0] for pos in trail1], [pos[1] for pos in trail1])
        line1.set_3d_properties([pos[2] for pos in trail1])

        # Linie 2
        if len(trail2) > MAX_LENGTH:
            trail2 = trail2[-MAX_LENGTH:]  # Behalte nur die letzten MAX_LENGTH Punkte
        line2.set_data([pos[0] for pos in trail2], [pos[1] for pos in trail2])
        line2.set_3d_properties([pos[2] for pos in trail2])

        # Linie 3
        if len(trail3) > MAX_LENGTH:
            trail3 = trail3[-MAX_LENGTH:]  # Behalte nur die letzten MAX_LENGTH Punkte
        line3.set_data([pos[0] for pos in trail3], [pos[1] for pos in trail3])
        line3.set_3d_properties([pos[2] for pos in trail3])

        # Optionale Planet-Linie
        if args.planet:
            if len(trail_planet) > MAX_LENGTH:
                trail_planet = trail_planet[-MAX_LENGTH:]  # Behalte nur die letzten MAX_LENGTH Punkte
            line_planet.set_data([pos[0] for pos in trail_planet], [pos[1] for pos in trail_planet])
            line_planet.set_3d_properties([pos[2] for pos in trail_planet])

    # Update der Punkte
    point1.set_data([x1], [y1])
    point1.set_3d_properties([z1])
    point2.set_data([x2], [y2])
    point2.set_3d_properties([z2])
    point3.set_data([x3], [y3])
    point3.set_3d_properties([z3])
    if args.planet:
        point_planet.set_data([x_planet], [y_planet])
        point_planet.set_3d_properties([z_planet])

    return (line1, line2, line3, line_planet, point1, point2, point3, point_planet) if args.planet and args.line else (point1, point2, point3, point_planet)

interval = 100 if args.slow else 20
ani = FuncAnimation(fig, update, init_func=init, frames=np.arange(0, 100000), interval=interval, blit=False)
plt.legend()
plt.show()
