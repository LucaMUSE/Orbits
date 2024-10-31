import numpy as np
import matplotlib.pyplot as plt

# Configuración de las condiciones iniciales
x0, y0, vx0, vy0 = 1.0, 0.0, 0.0, 1.0

# Intervalo de tiempo
t0, tf, dt = 0.0, 10.0, 0.01
N = int((tf - t0) / dt)  # Número de pasos de tiempo

# Creación del vector de tiempo
t_values = np.linspace(t0, tf, N)
n = len(t_values)

# Inicialización de las condiciones
U0 = np.array([x0, y0, vx0, vy0])
U_Euler = np.zeros((n, 4))
U_Euler[0] = U0
F_Euler = np.zeros(4)

# Método de EULER EXPLÍCITO
for i in range(n - 1):
    x, y, vx, vy = U_Euler[i]  # Extracción de las coordenadas

    # Cálculo de la fuerza
    norm_r = np.sqrt(x**2 + y**2)
    F_Euler[0], F_Euler[1] = vx, vy
    F_Euler[2] = -x / norm_r**3 if norm_r > 1e-10 else 0.0
    F_Euler[3] = -y / norm_r**3 if norm_r > 1e-10 else 0.0

    # Actualización de las posiciones y velocidades
    U_Euler[i + 1] = U_Euler[i] + dt * F_Euler

# Método de RUNGE-KUTTA de 4° orden
U_RK4 = np.zeros((n, 4))
U_RK4[0] = U0

for i in range(n - 1):
    x, y, vx, vy = U_RK4[i]
    norm_r = np.sqrt(x**2 + y**2)

    ax = -x / norm_r**3 if norm_r > 1e-10 else 0.0
    ay = -y / norm_r**3 if norm_r > 1e-10 else 0.0

    # Cálculo de los coeficientes k
    k1 = np.array([vx, vy, ax, ay])
    k2 = np.array([vx + 0.5 * dt * k1[2], vy + 0.5 * dt * k1[3], ax + 0.5 * dt * k1[2], ay + 0.5 * dt * k1[3]])
    k3 = np.array([vx + 0.5 * dt * k2[2], vy + 0.5 * dt * k2[3], ax + 0.5 * dt * k2[2], ay + 0.5 * dt * k2[3]])
    k4 = np.array([vx + dt * k3[2], vy + dt * k3[3], ax + dt * k3[2], ay + dt * k3[3]])

    U_RK4[i + 1] = U_RK4[i] + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

# Método de CRANK-NICKOLSON
U_CN = np.zeros((n, 4))
U_CN[0] = U0

for i in range(n - 1):
    x, y, vx, vy = U_CN[i]
    norm_r = np.sqrt(x**2 + y**2)

    ax = -x / norm_r**3 if norm_r > 1e-10 else 0.0
    ay = -y / norm_r**3 if norm_r > 1e-10 else 0.0

    # Actualización de las posiciones
    U_CN[i + 1, 0] = U_CN[i, 0] + dt * vx + 0.5 * dt**2 * ax
    U_CN[i + 1, 1] = U_CN[i, 1] + dt * vy + 0.5 * dt**2 * ay

    # Nuevas aceleraciones basadas en las posiciones actualizadas
    x_new, y_new = U_CN[i + 1, :2]
    norm_r_new = np.sqrt(x_new**2 + y_new**2)
    ax_new = -x_new / norm_r_new**3 if norm_r_new > 1e-10 else 0.0
    ay_new = -y_new / norm_r_new**3 if norm_r_new > 1e-10 else 0.0

    # Actualización de las velocidades
    U_CN[i + 1, 2] = U_CN[i, 2] + 0.5 * dt * (ax + ax_new)
    U_CN[i + 1, 3] = U_CN[i, 3] + 0.5 * dt * (ay + ay_new)

# Creación de las gráficas
plt.figure(figsize=(10, 6))
plt.plot(U_Euler[:, 0], U_Euler[:, 1], label="Euler Explícito", alpha=0.6)
plt.plot(U_RK4[:, 0], U_RK4[:, 1], label="Runge-Kutta 4 Orden", alpha=0.6)
plt.plot(U_CN[:, 0], U_CN[:, 1], label="Crank-Nickolson", alpha=0.6)
plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Solución de la EDO usando diferentes métodos numéricos")
plt.grid()
plt.show()
