import numpy as np

# 1. Euler Method
def Euler(U, t, h, F):
    return U + h * F(U, t)

# 2. Crank-Nicolson Method
def Crank_Nicolson(U, t, h, F, tol=1e-6, max_iter=10):
    U_new = U
    for _ in range(max_iter):
        U_new_next = U + h/2 * (F(U, t) + F(U_new, t + h))
        if np.linalg.norm(U_new_next - U_new) < tol:
            return U_new_next
        U_new = U_new_next
    return U_new

# 3. Runge-Kutta 4th Order (RK4) Method
def RK4(U, t, h, F):
    k1 = F(U, t)
    k2 = F(U + h/2 * k1, t + h/2)
    k3 = F(U + h/2 * k2, t + h/2)
    k4 = F(U + h * k3, t + h)
    return U + h/6 * (k1 + 2*k2 + 2*k3 + k4)

# 4. Inverse Euler Method
def Inverse_Euler(U, t, h, F, tol=1e-6, max_iter=10):
    U_new = U
    for _ in range(max_iter):
        U_new_next = U + h * F(U_new, t + h)
        if np.linalg.norm(U_new_next - U_new) < tol:
            return U_new_next
        U_new = U_new_next
    return U_new

# 5. Function to integrate a Cauchy problem
def integrate_Cauchy(U0, t0, t_end, h, F, method):
    U = U0
    t = t0
    while t < t_end:
        if method == 'Euler':
            U = Euler(U, t, h, F)
        elif method == 'Crank_Nicolson':
            U = Crank_Nicolson(U, t, h, F)
        elif method == 'RK4':
            U = RK4(U, t, h, F)
        elif method == 'Inverse_Euler':
            U = Inverse_Euler(U, t, h, F)
        t += h
    return U

# 6. Function to express the force of the Kepler movement
def F_Kepler(U, t):
    r = U[:2]
    r_dot = U[2:]
    r_norm = np.linalg.norm(r)
    r_double_dot = -r / r_norm**3
    return np.concatenate((r_dot, r_double_dot))

# 7. Integrate the Kepler problem
def integrate_Kepler(U0, t0, t_end, h, method):
    return integrate_Cauchy(U0, t0, t_end, h, F_Kepler, method)

# Example initial conditions and parameters
U0 = np.array([1, 0, 0, 1])  # Initial position and velocity
t0 = 0
t_end = 10
h = 0.01
methods = ['Euler', 'Crank_Nicolson', 'RK4', 'Inverse_Euler']

# Run integrations with different methods
results = {}
for method in methods:
    results[method] = integrate_Kepler(U0, t0, t_end, h, method)

# Print results
for method in results:
    print(f"Method: {method}, Final State: {results[method]}")

# 8. Increase and decrease the time step and explain the results
time_steps = [0.01, 0.05, 0.1]
for h in time_steps:
    print(f"Time step: {h}")
    for method in methods:
        result = integrate_Kepler(U0, t0, t_end, h, method)
        print(f"Method: {method}, Final State: {result}")
