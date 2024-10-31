from numpy import array, zeros, linspace  #array--vectores, linspace--partición equiespaciada
import matplotlib.pyplot as plt


#################################################### DATOS ########################################################
# Condiciones iniciales
x0 = 1.0
y0 = 0.0
vx0 = 0.0
vy0 = 1.0

# Instantes inicial, final y paso del tiempo
t0 = 0.0
tf = 10.0
dt = 0.01

# Número de intervalos (= nº de instantes de tiempo - 1)
N = int((tf - t0) / dt)


################################################### CÓDIGO ########################################################
# Inicializamos lado derecho de EDO para Euler 
# Vector instantes temporales
# Inicializamos el vector de intantes temporales en los que se calcula la solución
t_values = linspace(t0, tf, N)
n = len(t_values)

# Vector de condiciones iniciales
U0 = array ([x0, y0, vx0, vy0])

# Matrices de soluciones
U_Euler = zeros((n,4))
U_Euler[0] = U0 
F_Euler = zeros(4)


# Cálculos del EULER EXPLÍCITO
for i in range(n-1):
    # Coordenadas de posición y velocidad en cada intante de tiempo
    x = U_Euler[i, 0]
    y = U_Euler[i, 1]
    vx = U_Euler[i, 2]
    vy = U_Euler[i, 3]

    # Cálculo de la F_Euler
    F_Euler[0] = vx  # Componente de la velocidad en x
    F_Euler[1] = vy  # Componente de la velocidad en y

    # Cálculo de las otras dos componentes de F_Euler suponiendo 0 si se divide entre 0
    norm_r = (x**2 + y**2)**0.5
    if norm_r > 1e-10:
        F_Euler[2] = -x/norm_r**3
        F_Euler[3] = -y/norm_r**3
    else:
        F_Euler[2] = 0.0
        F_Euler[3] = 0.0
    
    # Cálculo de la U_Euler
    U_Euler[i+1, 0] = U_Euler[i][0] + dt * F_Euler[0]
    U_Euler[i+1, 1] = U_Euler[i][1] + dt * F_Euler[1]
    U_Euler[i+1, 2] = U_Euler[i][2] + dt * F_Euler[2]
    U_Euler[i+1, 3] = U_Euler[i][3] + dt * F_Euler[3]


# Cálculos del RUNGE-KUTTA DE 4 ETAPAS
U_RK4 = zeros((n, 4))
U_RK4[0] = U0

for i in range(n-1):
    x = U_RK4[i, 0]
    y = U_RK4[i, 1]
    vx = U_RK4[i, 2]
    vy = U_RK4[i, 3]

    # Cálculo de la aceleración evitando dividir entre 0
    norm_r = (x**2 + y**2)**0.5
    if norm_r > 1e-10:
        ax = -x/norm_r**3 
        ay = -y/norm_r**3  
    else:
        ax = 0.0
        ay = 0.0

    # Calculo de k1
    k1_x = vx
    k1_y = vy
    k1_vx = ax
    k1_vy = ay

    # Calculo de k2
    k2_x = vx + 0.5 * dt * k1_vx
    k2_y = vy + 0.5 * dt * k1_vy
    k2_vx = ax + 0.5 * dt * k1_vx
    k2_vy = ay + 0.5 * dt * k1_vy

    # Calculo de k3
    k3_x = vx + 0.5 * dt * k2_vx
    k3_y = vy + 0.5 * dt * k2_vy
    k3_vx = ax + 0.5 * dt * k2_vx
    k3_vy = ay + 0.5 * dt * k2_vy

    # Calculo de k4
    k4_x = vx + dt * k3_vx
    k4_y = vy + dt * k3_vy
    k4_vx = ax + dt * k3_vx
    k4_vy = ay + dt * k3_vy

    # Actualizar las posiciones y velocidades
    U_RK4[i+1, 0] = U_RK4[i, 0] + (dt/6) * (k1_x + 2*k2_x + 2*k3_x + k4_x)
    U_RK4[i+1, 1] = U_RK4[i, 1] + (dt/6) * (k1_y + 2*k2_y + 2*k3_y + k4_y)
    U_RK4[i+1, 2] = U_RK4[i, 2] + (dt/6) * (k1_vx + 2*k2_vx + 2*k3_vx + k4_vx)
    U_RK4[i+1, 3] = U_RK4[i, 3] + (dt/6) * (k1_vy + 2*k2_vy + 2*k3_vy + k4_vy)


# Cálculos del CRANK-NICKOLSON
U_CN = zeros((n, 4))
U_CN[0] = U0

for i in range(n-1):
    x = U_CN[i, 0]
    y = U_CN[i, 1]
    vx = U_CN[i, 2]
    vy = U_CN[i, 3]

    # Cálculo de la aceleración evitando dividir entre 0
    norm_r = (x**2 + y**2)**0.5
    if norm_r > 1e-10:
        ax = -x/norm_r**3 
        ay = -y/norm_r**3
    else:
        ax = 0.0
        ay = 0.0

    # Cálculo de las nuevas posiciones y velocidades
    U_CN[i+1, 0] = U_CN[i, 0] + dt * vx + 0.5 * dt**2 * ax
    U_CN[i+1, 1] = U_CN[i, 1] + dt * vy + 0.5 * dt**2 * ay
    
    # Aceleración en el nuevo tiempo (utilizando las posiciones actualizadas)
    x_new = U_CN[i+1, 0]
    y_new = U_CN[i+1, 1]
    norm_r_new = (x_new**2 + y_new**2)**0.5
    
    if norm_r_new > 1e-10:
        ax_new = -x_new / norm_r_new**3
        ay_new = -y_new / norm_r_new**3
    else:
        ax_new = 0.0
        ay_new = 0.0
    
    # Actualizar las velocidades usando la aceleración promedio
    U_CN[i+1, 2] = U_CN[i, 2] + 0.5 * dt * (ax + ax_new)
    U_CN[i+1, 3] = U_CN[i, 3] + 0.5 * dt * (ay + ay_new)


################################################# GRÁFICAS #########################################################
# Gráficas de las soluciones
plt.figure(figsize=(10, 6))
plt.plot(U_Euler[:, 0], U_Euler[:, 1], label="Euler Explícito", alpha=0.6)
plt.plot(U_RK4[:, 0], U_RK4[:, 1], label="Runge-Kutta 4 etapas", alpha=0.6)
plt.plot(U_CN[:, 0], U_CN[:, 1], label="Crank-Nickolson", alpha=0.6)
plt.legend()
plt.xlabel("x")
plt.ylabel("y")
plt.title("Solución de la EDO usando diferentes métodos numéricos")
plt.grid()
plt.show()