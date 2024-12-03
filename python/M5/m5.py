# Importar librerías necesarias
from numpy import array, reshape, zeros, linspace, pi, cos, sin
from numpy.linalg import norm
import matplotlib.pyplot as plt

################################################# FUNCIONES #######################################################
# Funcion para definir y guardar las condiciones iniciales en el vector U0
def Initial_conditions(Nb, Nc):
    '''
    INPUTS:
        - Nb: número de cuerpos
        - Nc: número de coordenadas
    '''
    U0 = zeros(2 * Nc * Nb) # Multiplicado por 2 porque cada cuerpo tiene 2 componentes
    U1 = reshape(U0, (Nb, Nc, 2)) # Matriz NbxNc y en cada coordenada un "vector" de 2 coordenadas
    r0 = reshape(U1[:, :, 0], (Nb, Nc))
    v0 = reshape(U1[:, :, 1], (Nb, Nc))

    # Parámetros para la configuración hexagonal (para que quede bonito)
    radius = 1.5  # Radio del hexágono
    angular_velocity = 0.3  # Velocidad angular para un movimiento circular

    for i in range(Nb):
        angle = 2 * pi * i / Nb  # Ángulo para distribuir los cuerpos
        r0[i, :] = [radius * cos(angle), radius * sin(angle), 0]
        v0[i, :] = [-radius * sin(angle) * angular_velocity, 
                     radius * cos(angle) * angular_velocity, 0]
    return U0 

# Función para el "N body problem"
def N_Body_Problem(U, t, Nb, Nc):
    '''
    INPUTS:
        - U: vector de condiciones iniciales
        - t: tiempo
        - Nb: número de cuerpos
        - Nc: número de coordenadas
    '''
    
    pu = reshape(U, (Nb, Nc, 2)) # Matriz NbxNc y en cada coordenada un "vector" de 2 coordenadas
    r = reshape(pu[:, :, 0], (Nb, Nc)) # Posiciones en pu
    v = reshape(pu[:, :, 1], (Nb, Nc)) # Velocidades en pu
    
    Fs = zeros(len(U))
    pFs = reshape(Fs, (Nb, Nc, 2))
    
    drdt = reshape(pFs[:, :, 0], (Nb, Nc)) # Derivadas de r en pFs
    dvdt = reshape(pFs[:, :, 1], (Nb, Nc)) # Derivadas de v en pFs
    
    dvdt[:,:] = 0 # Al inicio las aceleraciones son nulas (no las posiciones, drdt)
    
    for i in range(Nb): # Iterar sobre cada cuerpo
        drdt[i,:] = v[i,:]  # Cálculo de la velocidad de cada cuerpo
        for j in range(Nb):
            if j != i:  # Para que no se compare un cuerpo con él mismo
                d = r[j,:] - r[i,:]
                dvdt[i,:] = dvdt[i,:] +  d[:] / norm(d)**3
    return Fs

# Función para el problema de N cuerpos
def F_N_Body_Problem(U, t): 
    '''
    Esta función es como un "adaptador" para la función N_Body_Problem
    N_Body_Problem necesita 4 imputs, pero RK4 solo necesita 2, esta función adapta de uno a otro
    Sirve para meter los 2 imputs Nc, Nb necesarios
    INPUTS:
        - U: vector de condiciones iniciales
        - t: tiempo
        
    OUTPUT:
        - N_Body_Problem(U, t, Nb, Nc): función para el problema de N cuerpos
    '''
    
    return N_Body_Problem(U, t, Nb, Nc)

# Función para resolver el problema de Cauchy usando un esquema numérico
def Cauchy_problem(F, t, U0, Esquema):  
    N = len(t) - 1
    Nv = len(U0)
    U = zeros((N + 1, Nv))
    U[0, :] = U0
    
    for n in range(N): 
        U[n + 1, :] = Esquema(U[n, :], t[n + 1] - t[n], t[n], F) 
    return U

# Esquema de Runge-Kutta de 4 etapas
def RK4(U, dt, t, F):
    k1 = F(U, t)
    k2 = F(U + k1 * dt / 2, t + dt / 2)
    k3 = F(U + k2 * dt / 2, t + dt / 2)
    k4 = F(U + k3 * dt, t + dt)
    return U + (dt / 6) * (k1 + 2 * k2 + 2 * k3 + k4)

################################################### CÓDIGO ########################################################
# Definicion de los parametros de integracion
N = 1000 # Numero de pasos

# Definicion del numero de cuerpos y del numero de coordenadas por cuerpo
Nb = 3 # Número de cuerpos
Nc = 3 # Número de coordenadas

# Definicion del intervalo de tiempo
t0 = 0
tf = 4 * pi
t = linspace(t0, tf, N + 1) # Partición temporal en N+1 puntos

# Resolucion del prolema de Cauchy
U0 = Initial_conditions(Nb, Nc)
U = Cauchy_problem(F_N_Body_Problem, t, U0, RK4)

# Representacion de resultados
Us  = reshape(U, (N + 1, Nb, Nc, 2)) # Us[tiempo][cuerpo][coordenada][0, 1]
r   = reshape(Us[:, :, :, 0], (N + 1, Nb, Nc)) 

################################################# GRÁFICAS #########################################################
# Representacion de resultados (en 2 dimensiones)
for i in range(Nb):
    plt.plot(r[:, i, 0], r[:, i, 1])
plt.axis('equal')
plt.grid()
plt.show()

# Representacion de resultados (en 3 dimensiones)
Us = reshape(U, (Nb, Nc, 2, N + 1))
r = reshape(Us[:, :, 0, :], (Nb, Nc, N + 1))
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
for i in range(Nb):
    ax.plot3D(r[i, 0, :], r[i, 1, :], r[i, 2, :], label=f'Cuerpo {i + 1}')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.legend()
plt.show()
