"Solver for the 1D heat equation using the BTCS method."
"i'm using backward time, centered space (BTCS) method to solve the 1D heat equation."

import numpy as np

def initial_condition(x):
    """Initial condition function."""
    return x**2

def source_term(x, t):
    """Source term function."""
    return 0

def main():

    # Parameters and constants
    alpha = 1
    L = 1
    t_final=1
    n_nodos_t=3
    n_nodos_x=3

    delta_x = L/(n_nodos_x - 1)
    delta_t = L/(n_nodos_x - 1)

    # Fourier number
    F = alpha * delta_t / (delta_x ** 2)
    print("F =", F)

    x = np.linspace(0, L, n_nodos_x)
    t = np.linspace(0, t_final, n_nodos_t)

    T = np.zeros((n_nodos_t, n_nodos_x))
    T[0, :] = initial_condition(x)

    # Dirichlet boundary conditions

    T_0 = 0
    T_L = 1

    #
    A = np.zeros((n_nodos_x, n_nodos_x))

    for i in range(1, n_nodos_x - 1):
        A[i, i - 1] = -F
        A[i, i] = 1 + 2 * F
        A[i, i + 1] = -F

    A[0, 0] = 1
    A[-1, -1] = 1

    b = np.zeros(n_nodos_x)
    for n in range(0, n_nodos_t):
        for i in range(1, n_nodos_x - 1):
            b[i] = T[n-1, i] + delta_t + source_term(x[i], t[n])
        b[0] = T_0
        b[-1] = T_L

        T[n, :] = np.linalg.solve(A, b)
        print(f"t = {t[n]:.2f}, T = {T[n, :]}")


    pass