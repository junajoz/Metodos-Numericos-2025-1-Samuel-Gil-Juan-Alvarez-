"Solver for the Laplace equation for the electric potential in 2D using the finite difference method."

"""
Laplace 2D con FEM (P1 triangulares) en Ω=[0,8]x[0,6].
Dirichlet: V=0 en ∂Ω. Placas internas:
  y=2 (x∈[2,6]) → Vp1
  y=4 (x∈[2,6]) → Vp2
Autor: Juan Alvarez and Samuel Gil
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt

def solve_with_dirichlet_reduction(K, b, fixed_nodes, fixed_values):
    """
    Impone Dirichlet sin editar K:
      K_ff u_f = b_f - K_fF u_F
    Devuelve el vector solución completo u (incluye nodos fijos).
    """
    fixed_nodes  = np.asarray(fixed_nodes, dtype=int)
    fixed_values = np.asarray(fixed_values, dtype=float)
    N = K.shape[0]
    u = np.zeros(N)

    # DOFs libres
    all_idx = np.arange(N)
    free = np.setdiff1d(all_idx, fixed_nodes, assume_unique=False)

    # Sistema reducido
    K_ff = K[free][:, free]
    K_fF = K[free][:, fixed_nodes]
    b_f  = b[free] - K_fF @ fixed_values

    # Resolver y recomponer
    u[fixed_nodes] = fixed_values
    u[free] = spla.spsolve(K_ff, b_f)
    return u


def mesh_rectangular(Lx, Ly, nx, ny):
    """Genera una malla rectangular de nodos y elementos triangulares."""
    x = np.linspace(0, Lx, nx)
    y = np.linspace(0, Ly, ny)
    X, Y = np.meshgrid(x, y)
    nodes = np.column_stack((X.ravel(), Y.ravel()))
    rects = []
    for j in range(ny - 1):
        for i in range(nx - 1):
            n1 = j * nx + i
            n2 = n1 + 1
            n3 = n1 + nx
            n4 = n3 + 1
            rects.append([n1, n2, n4, n3])
    return x, y, X, Y, nodes, rects, nx, ny

def local_stiffness(xe, ye):
    """Matriz de rigidez local para elemento triangular P1."""
    Ae = 0.5 * np.linalg.det(np.array([[1, xe[0], ye[0]],
                                       [1, xe[1], ye[1]],
                                       [1, xe[2], ye[2]]]))
    b = np.array([ye[1] - ye[2], ye[2] - ye[0], ye[0] - ye[1]])
    c = np.array([xe[2] - xe[1], xe[0] - xe[2], xe[1] - xe[0]])
    ke = (np.outer(b, b) + np.outer(c, c)) / (4 * Ae)
    return ke, Ae

def assemble_global_matrix(nodes, elements):
    """Ensamblaje de la matriz global de rigidez."""
    n_nodes = nodes.shape[0]
    K = sp.lil_matrix((n_nodes, n_nodes))

    for elem in elements:
        xe = nodes[elem, 0]
        ye = nodes[elem, 1]
        ke, Ae = local_stiffness(xe, ye)
        for i in range(3):
            for j in range(3):
                K[elem[i], elem[j]] += ke[i, j]
    return K.tocsr()

def main():
    # Parámetros de la malla y del problema
    Lx, Ly = 8.0, 6.0
    nx, ny = 81, 61
    Vp1, Vp2 = 1.0, -1.0  # Potenciales en las placas internas
    x, y, X, Y, nodes, rects, nx, ny = mesh_rectangular(Lx, Ly, nx, ny)
    elements = []
    for r in rects:
        elements.append([r[0], r[1], r[2]])
        elements.append([r[0], r[2], r[3]])
    elements = np.array(elements)
    n_nodes = nodes.shape[0]
    K = assemble_global_matrix(nodes, elements)
    f = np.zeros(n_nodes)
    boundary_nodes = []
    boundary_values = []
    for i in range(n_nodes):
        xi, yi = nodes[i]
        if (np.isclose(xi, 0.0) or np.isclose(xi, Lx) or np.isclose(yi, 0.0) or np.isclose(yi, Ly)):
            boundary_nodes.append(i)
            boundary_values.append(0.0)
        elif 2.0 <= xi <= 6.0 and np.isclose(yi, 2.0):
            boundary_nodes.append(i)
            boundary_values.append(Vp1)
        elif 2.0 <= xi <= 6.0 and np.isclose(yi, 4.0):
            boundary_nodes.append(i)
            boundary_values.append(Vp2)
    V = solve_with_dirichlet_reduction(K, f, boundary_nodes, boundary_values)

    V_grid = V.reshape((ny, nx))
    plt.figure(figsize=(8, 6))
    cp = plt.contourf(X, Y, V_grid, levels=50, cmap='jet')
    plt.colorbar(cp)
    plt.title('Potencial eléctrico V en 2D')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    plt.show()
    # Guardar resultados en archivo
    np.savetxt('laplace2d_solution.txt', np.column_stack((nodes, V)), header='x y V', comments='')

if __name__ == "__main__":
    main()