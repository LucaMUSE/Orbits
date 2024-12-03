from numpy import array, zeros, linspace, concatenate
from numpy.linalg import norm
import matplotlib.pyplot as plt

################################################## FUNCIONES #######################################################

# KEPLER
def Kepler(U): #U = vector
    r = U[0:2]
    rdot = U[2:4]
    F = concatenate ( (rdot, -r/norm(r)**3), axis =0)
    
    return F


# Esquema EULER EXPLÍCITO
def Euler (F, U, dt, t):
    
    return U + dt * F(U,t)


def Cauchy(Esquema, F, U0, t): 
    N = len(t)-1  # Por como se empieza a contar en Python
    U = zeros((N+1,len(U0)))
    
    U[0,:] = U0
    for n in range(0,N):
        U[n+1,:] = Esquema ( Euler, U[n,:], t[n+1]-t[n], t[n] )
    return

# Esquema RUNGE-KUTTA órden 4
def RK4 (F, U, dt, t):
    k1 = F(U, t)
    k2 = F ( U + k1 * dt/2, t + dt/2)
    k3 = F ( U + k2 * dt/2, t + dt/2)
    k4 = F ( U + k3 * dt , t + dt/2)
    return

# EULER para un oscilador armónico
def Oscilador (U, t):
    xdot = U[0]
    x = U[1]
    F = array ([xdot, -x])
    return




################################################### CÓDIGO #########################################################
F_Euler = zeros[n,n]
U_Euler = zeros[n,n]

F_Euler[n,:] = Kepler(U_Euler[n,:])

U_Euler[n+1,:] = U_Euler[n,:] + (t[n+1]-t[n]) * F_Euler[n,:]

U_Euler[n+1,:] = Euler (Kepler, U_Euler[n,:], dt, t[n])  

# CAUCHY
U_Euler = Cauchy (Euler, Kepler, U0, t)

# RK4
U_Euler = Cauchy (Euler, Kepler, U0, t)

# oscilador armónico
U0_osc = ([x0_osc, vx0_osc])
x0_osc = 1
vx0_osc = 0
U_osc_Euler = Cauchy (Euler, Oscilador, U0_osc, t)