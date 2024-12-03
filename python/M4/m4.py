# Importar librerías necesarias
from numpy import zeros, linspace, abs, float64
import matplotlib.pyplot as plt

################################################# FUNCIONES #######################################################
# REGIÓN DE ESTABILIDAD GENERAL
def Region_Estabilidad(Scheme, N, x0, xf, y0, yf):
    '''
    INPUTS:
        - Scheme: esquema numérico a utilizar
        - N: número de puntos en la malla
        - x0, xf: intervalo en x
        - y0, yf: intervalo en y
    '''
    
    x = linspace(x0, xf, N)  # Malla en x
    y = linspace(y0, yf, N)  # Malla en y
    rho = zeros((N, N), dtype=float64)  # Matriz de estabilidad
    F = lambda u, t: w * u  # Función anónima que representa la ecuación diferencial a resolver
    
    for i in range(N):
        for j in range(N):
            w = complex(x[i], y[j])  # Número complejo
            r = Scheme(1.0, 1.0, 0.0, F)  # Paso temporal, valor inicial, tiempo inicial, función a resolver
            rho[i, j] = abs(r)
    return rho, x, y

# Esquemas numéricos
def Euler(U, dt, t, F):
    return U + dt * F(U, t)

def RK4(U, dt, t, F):
    k1 = F(U, t)
    k2 = F(U + k1 * dt / 2, t + dt / 2)
    k3 = F(U + k2 * dt / 2, t + dt / 2)
    k4 = F(U + k3 * dt, t + dt)
    return U + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

def Crank_Nickolson(U, dt, t, F):
    def G(X):  # Siendo X = U(n-1)
        return X - U - dt / 2 * (F(X, t) + F(U, t))
    from scipy.optimize import newton
    return newton(G, U)

def Euler_Inverso(U, dt, t, F):
    def G(X):  # Siendo X = U(n-1)
        return X - U - dt * F(X, t)
    from scipy.optimize import newton
    return newton(G, U)

# REGIÓN DE ESTABILIDAD DE EULER
def Region_Estabilidad_Euler(N, x0, xf, y0, yf):
    w = linspace(x0, xf, N) + 1j * linspace(y0, yf, N)[:, None]
    r = 1 + w
    return abs(r) <= 1

# REGIÓN DE ESTABILIDAD DE RK4
def Region_Estabilidad_RK4(N, x0, xf, y0, yf):
    w = linspace(x0, xf, N) + 1j * linspace(y0, yf, N)[:, None]
    r = 1 + w + w**2 / 2 + w**3 / 6 + w**4 / 24
    return abs(r) <= 1

# REGIÓN DE ESTABILIDAD DE CRANK-NICKOLSON
def Region_Estabilidad_Crank_Nickolson(N, x0, xf, y0, yf):
    w = linspace(x0, xf, N) + 1j * linspace(y0, yf, N)[:, None]
    r = (1 + w / 2) / (1 - w / 2)
    return abs(r) <= 1

# REGIÓN DE ESTABILIDAD DE EULER INVERSO
def Region_Estabilidad_Euler_Inverso(N, x0, xf, y0, yf):
    w = linspace(x0, xf, N) + 1j * linspace(y0, yf, N)[:, None]
    r = 1 / (1 - w)
    return abs(r) <= 1

############################################ CÓDIGO & GRÁFICAS ####################################################
# Graficar las regiones de estabilidad
plt.figure(figsize=(18, 6))

# Región de estabilidad de Euler
plt.subplot(1, 4, 1)
plt.imshow(Region_Estabilidad_Euler(400, -3, 3, -3, 3), extent=[-3, 3, -3, 3], origin='lower', cmap='Greys')
plt.title('Región de estabilidad de Euler')
plt.xlabel('Re(r)')
plt.ylabel('Im(r)')

# Región de estabilidad de RK4
plt.subplot(1, 4, 2)
plt.imshow(Region_Estabilidad_RK4(400, -3, 3, -3, 3), extent=[-3, 3, -3, 3], origin='lower', cmap='Greys')
plt.title('Región de estabilidad de RK4')
plt.xlabel('Re(r)')
plt.ylabel('Im(r)')

# Región de estabilidad de Crank-Nicolson
plt.subplot(1, 4, 3)
plt.imshow(Region_Estabilidad_Crank_Nickolson(400, -3, 3, -3, 3), extent=[-3, 3, -3, 3], origin='lower', cmap='Greys')
plt.title('Región de estabilidad de Crank-Nicolson')
plt.xlabel('Re(r)')
plt.ylabel('Im(r)')

# Región de estabilidad de Euler Inverso
plt.subplot(1, 4, 4)
plt.imshow(Region_Estabilidad_Euler_Inverso(400, -3, 3, -3, 3), extent=[-3, 3, -3, 3], origin='lower', cmap='Greys')
plt.title('Región de estabilidad de Euler Inverso')
plt.xlabel('Re(r)')
plt.ylabel('Im(r)')

plt.tight_layout()
plt.show()
