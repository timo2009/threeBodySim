import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.animation import FuncAnimation

# Gravitationskonstante
G = 1.0

# Massen der drei Körper
m1, m2, m3 = 1.0, 1.0, 1.0

# Anfangsbedingungen: Positionen und Geschwindigkeiten der drei Körper
def initial_conditions():
    # Positionen
    r1 = [-5.0, 0.0, 0.0]
    r2 = [5.0, 0.0, 0.0]
    r3 = [0.0, 5.0, 0.0]

    # Geschwindigkeiten
    v1 = [-0.1, 0.0, 0.1]
    v2 = [0.0, 0.1, -0.1]
    v3 = [0.1, -0.1, 0.0]

    # Kombinierte Anfangsbedingungen
    return np.array(r1 + r2 + r3 + v1 + v2 + v3)

# Funktion, die die Bewegungsgleichungen beschreibt
def derivatives(state, t, G, m1, m2, m3):
    # Entpacke den Zustand
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

    # Rückgabe der Ableitungen (dx/dt, dv/dt)
    return np.array([
        vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3,  # Geschwindigkeiten
        a1[0], a1[1], a1[2], a2[0], a2[1], a2[2], a3[0], a3[1], a3[2]  # Beschleunigungen
    ])

# Zeitspanne der Simulation
t = np.linspace(0, 5000, 50000)

# Anfangsbedingungen laden
initial_state = initial_conditions()

# Numerische Integration
solution = odeint(derivatives, initial_state, t, args=(G, m1, m2, m3))

# Ergebnisse extrahieren
x1, y1, z1 = solution[:, 0], solution[:, 1], solution[:, 2]
x2, y2, z2 = solution[:, 3], solution[:, 4], solution[:, 5]
x3, y3, z3 = solution[:, 6], solution[:, 7], solution[:, 8]

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
line1, = ax.plot([], [], [], lw=1, label="Körper 1", color="blue")
line2, = ax.plot([], [], [], lw=1, label="Körper 2", color="green")
line3, = ax.plot([], [], [], lw=1, label="Körper 3", color="orange")
point1, = ax.plot([], [], [], 'o', color="blue")
point2, = ax.plot([], [], [], 'o', color="green")
point3, = ax.plot([], [], [], 'o', color="orange")

# Speicher für die Bahnen
trail1, trail2, trail3 = [], [], []

# Initialisierungsfunktion
def init():
    line1.set_data([], [])
    line1.set_3d_properties([])
    line2.set_data([], [])
    line2.set_3d_properties([])
    line3.set_data([], [])
    line3.set_3d_properties([])
    point1.set_data([], [])
    point1.set_3d_properties([])
    point2.set_data([], [])
    point2.set_3d_properties([])
    point3.set_data([], [])
    point3.set_3d_properties([])
    return line1, line2, line3, point1, point2, point3

# Update-Funktion für die Animation
def update(frame):
    global trail1, trail2, trail3

    # Positionen der Körper
    trail1.append((x1[frame], y1[frame], z1[frame]))
    trail2.append((x2[frame], y2[frame], z2[frame]))
    trail3.append((x3[frame], y3[frame], z3[frame]))

    # Update der Linien und Punkte
    line1.set_data([pos[0] for pos in trail1], [pos[1] for pos in trail1])
    line1.set_3d_properties([pos[2] for pos in trail1])
    line2.set_data([pos[0] for pos in trail2], [pos[1] for pos in trail2])
    line2.set_3d_properties([pos[2] for pos in trail2])
    line3.set_data([pos[0] for pos in trail3], [pos[1] for pos in trail3])
    line3.set_3d_properties([pos[2] for pos in trail3])

    # Punktpositionen als Liste
    point1.set_data([x1[frame]], [y1[frame]])  # Einzeln als Liste
    point1.set_3d_properties([z1[frame]])  # Einzeln als Liste
    point2.set_data([x2[frame]], [y2[frame]])  # Einzeln als Liste
    point2.set_3d_properties([z2[frame]])  # Einzeln als Liste
    point3.set_data([x3[frame]], [y3[frame]])  # Einzeln als Liste
    point3.set_3d_properties([z3[frame]])  # Einzeln als Liste

    return line1, line2, line3, point1, point2, point3

# Erstellen der Animation
ani = FuncAnimation(fig, update, frames=len(t), init_func=init, blit=False, interval=20)  # interval erhöht für flüssige Animation
plt.legend()
plt.show()