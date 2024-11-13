import numpy as np
import matplotlib.pyplot as plt

# Condiciones iniciales para la órbita de Kepler
x_kepler, y_kepler, xd_kepler, yd_kepler = 1, 0, 0, 1
t0, tf, N = 0, 20, 200
deltat = (tf - t0) / N
U0_kepler = np.array([x_kepler, y_kepler, xd_kepler, yd_kepler])

# Función para calcular las derivadas del sistema (rhs) para la órbita de Kepler
def kepler_rhs(U):
    x, y, xd, yd = U
    r3 = (x**2 + y**2)**(3/2)  # Cálculo de la distancia al cubo
    return np.array([xd, yd, -x/r3, -y/r3])  # Ecuaciones del movimiento en el espacio

# Método de Euler explícito
def euler_step(U, deltat):
    return U + deltat * kepler_rhs(U)  # Avance en el tiempo según Euler explícito

# Método de Runge-Kutta 2 pasos
def rk2_step(U, deltat, t):
    k1 = kepler_rhs(U)  # Cálculo del primer paso de Runge-Kutta
    k2 = kepler_rhs(U + deltat * k1 / 2)  # Cálculo del segundo paso de Runge-Kutta
    return U + deltat * k2  # Avance en el tiempo usando Runge-Kutta 2

# Método de Runge-Kutta 3 pasos
def rk3_step(U, deltat, t):
    k1 = kepler_rhs(U)  # Primer paso de Runge-Kutta 3
    k2 = kepler_rhs(U + deltat * k1 / 3)  # Segundo paso
    k3 = kepler_rhs(U + deltat * k2 * (2 / 3))  # Tercer paso
    return U + deltat * (k1 + 3 * k3) / 4  # Avance en el tiempo con Runge-Kutta 3

# Método de Runge-Kutta 4 pasos
def rk4_step(U, deltat, t):
    k1 = kepler_rhs(U)  # Primer paso de Runge-Kutta 4
    k2 = kepler_rhs(U + deltat * k1 / 2)  # Segundo paso
    k3 = kepler_rhs(U + deltat * k2 / 2)  # Tercer paso
    k4 = kepler_rhs(U + deltat * k3)  # Cuarto paso
    return U + deltat * (k1 + 2*k2 + 2*k3 + k4) / 6  # Avance en el tiempo con Runge-Kutta 4

# Inicialización de los arrays para almacenar los resultados
U_euler = np.zeros((N+1, 4))
U_rk2 = np.zeros((N+1, 4))
U_rk3 = np.zeros((N+1, 4))
U_rk4 = np.zeros((N+1, 4))

U_euler[0, :] = U0_kepler  # Condición inicial para Euler
U_rk2[0, :] = U0_kepler  # Condición inicial para RK2
U_rk3[0, :] = U0_kepler  # Condición inicial para RK3
U_rk4[0, :] = U0_kepler  # Condición inicial para RK4

# Tiempo
t = np.linspace(t0, tf, N+1)

# Cálculo de las trayectorias usando diferentes métodos
for n in range(N):
    U_euler[n+1, :] = euler_step(U_euler[n, :], deltat)  # Avance de Euler
    U_rk2[n+1, :] = rk2_step(U_rk2[n, :], deltat, t[n])  # Avance de Runge-Kutta 2
    U_rk3[n+1, :] = rk3_step(U_rk3[n, :], deltat, t[n])  # Avance de Runge-Kutta 3
    U_rk4[n+1, :] = rk4_step(U_rk4[n, :], deltat, t[n])  # Avance de Runge-Kutta 4

# Gráficas de las órbitas
plt.figure(figsize=(8, 5))
plt.axis("equal")
plt.plot(U_euler[:, 0], U_euler[:, 1], '-b', lw=1, label="Euler explícito")
plt.plot(U_rk2[:, 0], U_rk2[:, 1], '-r', lw=1, label="Runge-Kutta 2")
plt.plot(U_rk3[:, 0], U_rk3[:, 1], '--g', lw=1, label="Runge-Kutta 3")
plt.plot(U_rk4[:, 0], U_rk4[:, 1], ':y', lw=1, label="Runge-Kutta 4")
plt.legend()
plt.xlabel('Coordenada x')
plt.ylabel('Coordenada y')
plt.title(r'Órbita de Kepler ($\Delta t$ = {})'.format(round(deltat, 2)))
plt.grid()
plt.show()
