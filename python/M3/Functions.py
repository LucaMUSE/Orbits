from numpy import concatenate, zeros, linspace, array, exp, log10
from numpy.linalg import norm
from scipy.optimize import newton



def Problema_Cauchy(Esquema, F, U0, t):
    '''''''''''
    --INPUTS--

    Esquema(U, F, t): función que representa el esquema numérico a utilizar
    F(U,t): función a resolver
    U0: vector de condiciones iniciales
    t: partición temporal
    '''''''''''



    N = len(t) - 1
    U = zeros([N+1, len(U0)])
    U[0,:] = U0

    for n in range(N):

        U[n+1,:] = Esquema(U[n,:], F, t)

    return U

def Kepler(U,t):
    '''''''''''
    --INPUTS--

    U: vector de estado (posición, velocidad)
    t: partición temporal 
    '''''''''''


    F = concatenate((U[2:4],-U[0:2]/norm(U[0:2])**3), axis = 0)

    return F

def Oscilador(U,t):
    '''''''''''
    --INPUTS--

    U: vector de estado (posición, velocidad)
    t: partición temporal 
    '''''''''''

    F = array([ U[1], -U[0]])

    return F




def Euler_Explicito(U, F, t):
    '''''''''''
    --INPUTS--

    U: vector de estado (posición, velocidad)
    F(U,t): función a resolver
    t: partición temporal 
    '''''''''''

    dt = t[1]-t[0]

    return U + dt*F(U,t)

def RK4(U, F, t):
    '''''''''''
    --INPUTS--

    U: vector de estado (posición, velocidad)
    F(U,t): función a resolver
    t: partición temporal 
    '''''''''''

    dt = t[1]-t[0]

    k1 = F(U,t)
    k2 = F(U + dt*k1/2, t + dt/2)
    k3 = F(U + dt*k2/2, t + dt/2)
    k4 = F(U + dt*k3, t + dt)

    return U + dt*(k1 + 2*k2 + 2*k3 + k4)/6

def Euler_Implicito(U, F, t):
    '''''''''''
    --INPUTS--

    U: vector de estado (posición, velocidad)
    F(U,t): función a resolver
    t: partición temporal 
    '''''''''''

    dt = t[1]-t[0]

    def G(X):

        return X - U - dt*F(X, t)
    
    return newton(G, U, maxiter = 150) #Nota: Escribir Newton hecho por mí.

def Crank_Nicolson(U, F, t):
    '''''''''''
    --INPUTS--

    U: vector de estado (posición, velocidad)
    F(U,t): función a resolver
    t: partición temporal
    '''''''''''

    dt = t[1]-t[0]

    def G(X):

        return X - U - dt*(F(X, t)+F(U,t))
    
    return newton(G, U)

def Problem_Error(U0, F, Problema, Esquema, t, q):  #Este método esta escrito sólo para usarse con dt1 = dt2/2, generalizar en el futuro
    '''''''''''
    --INPUTS--

    U0: vector de condiciones iniciales
    F(U,t): función a resolver
    Problema(Esquema, F, U0, t): Función que representa el problema a resolver (Cauchy hasta el momento)
    Esquema(U, F, t): función que representa el esquema numérico a utilizar
    t: partición temporal 
    '''''''''''

    N = len(t)-1
    t1 = t
    t2 = linspace(t[0], t[-1], 2*N+1) #Refinamiento de malla

    Error = zeros([len(t1), len(U0)])           

    U1 = Problema(Esquema, F, U0, t1) #Solución al problema con la malla original
    U2 = Problema(Esquema, F, U0, t2) #Solución al problema con la malla refinada

    for i in range(len(t)):
        Error[i,:] = (U2[2*i,:] - U1[i,:]) / (1 - 1/2**q)

    return Error

def Problem_Error_Convergencia(U0, F, Problema, Esquema, t):  

    N = len(t)-1
    t1 = t
    t2 = linspace(t[0], t[-1], 2*N+1) #Refinamiento de malla

    Error = zeros([len(t1), len(U0)])           

    U1 = Problema(Esquema, F, U0, t1) #Solución al problema con la malla original
    U2 = Problema(Esquema, F, U0, t2) #Solución al problema con la malla refinada

    for i in range(len(t)):
        Error[i,:] = (U2[2*i,:] - U1[i,:])

    return Error

def Convergencia(U0, F, Error, Problema, Esquema, t):
    '''''''''''
    --INPUTS--
    
        U0: Vector del estado inicial
        F: Función a resolver
        Error(U0, F, Problema, Esquema, t): Función que devuelve un vector con el error de un esquema en cada paso temporal
        Esquema: Esquema temporal a resolver
        t: partición temporal 

    '''''''''''
    np = 15 #Número de puntos de la regresión 
    logE = zeros(np)
    logN = zeros(np)
    N = len(t-1)
    t1 = t
    for i in range(np):
        
        E = Error(U0, F, Problema, Esquema, t1)
        logE[i] = log10(norm(E[-1,:]))
        logN[i] = log10(N)
        N = 2*N
        t1 = linspace(t[0], t[-1], N+1)    

    return logN, logE


def Newton(F, x_0, Fprima = None, tol = 1e-8, maxiter=50):
    '''''''''''
    
    INPUTS
        F: Función escalar de la que sacar las raíces
        x_0: Punto inicial del eje x en el que se comienza la iteración  
        Fprima: derivada de F. (Si no se introduce se calcula dentro de la función)
        tol: tolerancia (por defecto es 10e-8)
        maxiter: número máximo de iteraciones

    
    '''''''''''
    def Fp(x):
        if Fprima == None:
            delta = 1e-4
            return (F(x+delta)-F(x-delta))/(2*delta)
        else:
            return Fprima(x)
   
   
    xn = x_0
    Error = tol + 1
    iter = 0

    while Error > tol and iter < maxiter:

        xn1 = xn - F(xn)/Fp(xn)
        Error = abs(xn-xn1)
        xn = xn1

        iter += 1
        print('Error:', Error)
    print('Número de iteraciones: ', iter)
    return xn

def jorge(x):
    return exp(x)-2*x-2 # Al usar la exp de numpy, se aceptan inputs tanto como escalar o vectorial


def dif_jorge(x):
    return exp(x)-2 





