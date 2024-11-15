from Functions import Problema_Cauchy, Kepler, Oscilador, Euler_Explicito, RK4, Euler_Implicito, Crank_Nicolson, Problem_Error, Problem_Error_Convergencia, Convergencia
from numpy import array, linspace
import matplotlib.pyplot as plt
from scipy.stats import linregress

U0 = array([1, 0, 0, 1])
t0 = 0
tf = 20
N = 500
t = linspace(t0, tf, N+1)

input1 = input("Introduce el esquema numérico deseado(EE, RK4, EI, CN): ")
input2 = input("Introduce qué desea saber(Error, Conver)")

if input1 == "EE":
    esquema = Euler_Explicito
    q = 1

elif input1 == "RK4":
    esquema = RK4
    q = 4

elif input1 == "EI":
    esquema = Euler_Implicito
    q = 1

elif input1 == "CN":
    esquema = Crank_Nicolson
    q = 2
else:
    print("Esquema no válido")

if input2 == "Error":

    Error = Problem_Error(U0, Kepler, Problema_Cauchy, esquema, t, q)

    plt.axis('equal') 
    plt.xlabel('tiempo')
    plt.ylabel('Error')
    plt.plot(t, Error[:,0], '-b', label = 'X')
    plt.plot(t, Error[:,1], '-r', label = 'Y')
    plt.plot(t, Error[:,2], '--m', label = 'Vx')
    plt.plot(t, Error[:,3], '--c', label = 'Vy')
    plt.legend()
    plt.show()

elif input2 == "Conver":

    logN, logE = Convergencia(U0, Kepler, Problem_Error_Convergencia, Problema_Cauchy, esquema, t) #Ajustar N en función del esquema para ver únicamente la parte recta

    plt.axis('equal') 
    plt.xlabel('logN')
    plt.ylabel('logE')
    plt.plot(logN, logE, '-b')
    plt.show()

    # q = -(logE[-1]-logE[-2])/(logN[-1]-logN[-2])

    # print('El orden de', input1, 'es', q)
