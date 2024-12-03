# Importar librerías necesarias
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import newton

#################################################### DATOS ########################################################
# Condiciones iniciales
x0, y0 = 1.0, 0.0
vx0, vy0 = 0.0, 1.0

# Variables de tiempo y número de intervalos 
tf = 7
N = 2000

################################################# FUNCIONES #######################################################
# Función que define el problema de Kepler
def Kepler(U, t):
    x, y, dxdt, dydt = U
    d = (x**2 + y**2)**1.5
    return np.array([dxdt, dydt, -x/d, -y/d])

# Problema de Cauchy: obtener la solución de un problema de CI dada una CI y un esquema temporal
def Cauchy_problem(F, t, U0, Esquema):
    N = len(t) - 1
    Nv = len(U0)
    U = np.zeros((N + 1, Nv), dtype=float) 
    U[0, :] = U0
    
    for n in range(N):
        U[n + 1, :] = Esquema(U[n, :], t[n + 1] - t[n], t[n], F)
    return U

# Esquema de Euler Explícito
def Euler(U, dt, t, F):
    return U + dt * F(U, t)

# Esquema de Runge-Kutta de órden 4
def RK4(U, dt, t, F):
    k1 = F(U, t)
    k2 = F(U + k1 * dt / 2, t + dt / 2)
    k3 = F(U + k2 * dt / 2, t + dt / 2)
    k4 = F(U + k3 * dt, t + dt)
    return U + (dt / 6) * (k1 + 2*k2 + 2*k3 + k4)

# Esquema implícito Crank-Nickolson
def Crank_Nickolson(U, dt, t, F):
    def G(X):  # Siendo X = U(n+1)
        return X - U - dt / 2 * (F(X, t) + F(U, t))
    return newton(G, U)

# Esquema implícito Euler Inverso
def Euler_Inverso(U, dt, t, F):
    def G(X):  # Siendo X = U(n+1)
        return X - U - dt * F(X, t)
    return newton(G, U)

################################################### CÓDIGO ########################################################
# Separación equiespaciada de instantes de tiempo en los que calcular la solución
t = np.linspace(0, tf, N)

# Vector de condiciones iniciales
U0 = np.array([x0, y0, vx0, vy0])

# Resolución del problema de Cauchy con diferentes métodos
U_Euler = Cauchy_problem(Kepler, t, U0, Euler)
U_RK4 = Cauchy_problem(Kepler, t, U0, RK4)
U_CN = Cauchy_problem(Kepler, t, U0, Crank_Nickolson)
U_EulerInverso = Cauchy_problem(Kepler, t, U0, Euler_Inverso)

################################################# GRÁFICAS ########################################################
# Gráficas de las soluciones obtenidas
plt.figure(figsize=(10, 6))
plt.plot(U_Euler[:, 0], U_Euler[:, 1], label="Euler Explícito", alpha=0.6)
plt.plot(U_RK4[:, 0], U_RK4[:, 1], label="Runge-Kutta 4 etapas", alpha=0.6)
plt.plot(U_CN[:, 0], U_CN[:, 1], label="Crank-Nickolson", alpha=0.6)
plt.plot(U_EulerInverso[:, 0], U_EulerInverso[:, 1], label="Euler Inverso", alpha=0.6)
plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Solución de la EDO usando diferentes métodos numéricos")
plt.grid()
plt.show()
