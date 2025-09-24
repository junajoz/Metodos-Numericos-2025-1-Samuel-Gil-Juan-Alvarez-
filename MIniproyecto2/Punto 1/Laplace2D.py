"Laplace equation for the electric potential in 2D using the finite difference method."

"""
Laplace 2D con FEM (P1 triangulares) en Ω=[0,8]x[0,6].
Condiciones de Dirichlet:
  - V=0 en la frontera ∂Ω
  - Placa inferior:  y=2,  x∈[2,6] → Vp1
  - Placa superior:  y=4,  x∈[2,6] → Vp2

Autor: Juan Álvarez y Samuel Gil
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt


# ----------------------------- Utilidades FEM -----------------------------

def solve_with_dirichlet_reduction(K, b, fixed_nodes, fixed_values):
    """
    Impone Dirichlet SIN editar la matriz K:
      Sea F el conjunto de nodos fijos y f el conjunto libre.
      K_ff u_f = b_f - K_fF u_F
    Retorna el vector solución completo u (incluye nodos fijos).
    """
    fixed_nodes  = np.asarray(fixed_nodes, dtype=int)
    fixed_values = np.asarray(fixed_values, dtype=float)
    N = K.shape[0]
    u = np.zeros(N)

    all_idx = np.arange(N)
    free = np.setdiff1d(all_idx, fixed_nodes, assume_unique=False)

    K_ff = K[free][:, free]
    K_fF = K[free][:, fixed_nodes]
    b_f  = b[free] - K_fF @ fixed_values

    u[fixed_nodes] = fixed_values
    u[free] = spla.spsolve(K_ff, b_f)
    return u


def mesh_rectangular(Lx, Ly, nx, ny):
    """
    Genera una malla rectangular de nodos y la lista de rectángulos;
    luego, cada rectángulo se parte en dos triángulos para usar P1.
    """
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
            rects.append([n1, n2, n4, n3])   # [ab-izq, ab-der, ar-der, ar-izq]
    return x, y, X, Y, nodes, rects, nx, ny


def local_stiffness(xe, ye):
    """
    Matriz de rigidez local 3x3 para un triángulo lineal (P1).
    ke_ij = (b_i b_j + c_i c_j) / (4 A_e)
    """
    Ae = 0.5 * np.linalg.det(np.array([[1, xe[0], ye[0]],
                                       [1, xe[1], ye[1]],
                                       [1, xe[2], ye[2]]]))
    b = np.array([ye[1] - ye[2], ye[2] - ye[0], ye[0] - ye[1]])
    c = np.array([xe[2] - xe[1], xe[0] - xe[2], xe[1] - xe[0]])
    ke = (np.outer(b, b) + np.outer(c, c)) / (4.0 * Ae)
    return ke


def assemble_global_matrix(nodes, elements):
    """
    Ensamblaje de la matriz global de rigidez K (formato LIL→CSR).
    """
    n_nodes = nodes.shape[0]
    K = sp.lil_matrix((n_nodes, n_nodes))
    for elem in elements:
        xe = nodes[elem, 0]
        ye = nodes[elem, 1]
        ke = local_stiffness(xe, ye)
        for i in range(3):
            for j in range(3):
                K[elem[i], elem[j]] += ke[i, j]
    return K.tocsr()


# ------------------------ Visualización y guardado ------------------------

def trazar_potencial_y_campo(X, Y, V_grid, Vp1, Vp2,
                             ruta_salida="laplace2d.png",
                             dibujar_placas=True):
    """
    Dibuja equipotenciales + mapa de color del potencial y, además,
    líneas de campo eléctrico E = -∇V usando streamplot. Guarda PNG.
    """
    # Gradientes con espaciamientos reales (filas → y, columnas → x)
    dx = X[0, 1] - X[0, 0]
    dy = Y[1, 0] - Y[0, 0]
    dV_dy, dV_dx = np.gradient(V_grid, dy, dx)
    Ex = -dV_dx
    Ey = -dV_dy

    fig, ax = plt.subplots(figsize=(8, 6))
    cf = ax.contourf(X, Y, V_grid, levels=50, cmap="jet")
    ax.contour(X, Y, V_grid, levels=18, linewidths=0.3, colors="k", alpha=0.5)
    # líneas de campo
    sp = ax.streamplot(X, Y, Ex, Ey, density=1.2, linewidth=0.6,
                    color="k", arrowsize=0.6)
    sp.lines.set_alpha(0.5)   # líneas de flujo
    sp.arrows.set_alpha(0.5)  # puntas de flecha

    if dibujar_placas:
        ax.plot([2, 6], [2, 2], linewidth=3, color="w")   # placa Vp1
        ax.plot([2, 6], [4, 4], linewidth=3, color="w")   # placa Vp2

    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label("Potencial [V]")

    ax.set_title("Potencial y líneas de campo (Laplace 2D)")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    plt.tight_layout()
    plt.savefig(ruta_salida, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(fig)  # cierra la figura para no duplicar en notebooks


# --------------------------------- Main -----------------------------------

def main():
    # Parámetros del dominio y de la malla
    Lx, Ly = 8.0, 6.0
    nx, ny = 81, 61
    Vp1, Vp2 = 1.0, -1.0  # Potenciales en las placas

    # 1) Malla
    x, y, X, Y, nodes, rects, nx, ny = mesh_rectangular(Lx, Ly, nx, ny)

    # 2) Conectividad triangular (dos triángulos por rectángulo)
    elements = []
    for r in rects:
        elements.append([r[0], r[1], r[2]])
        elements.append([r[0], r[2], r[3]])
    elements = np.array(elements, dtype=int)

    # 3) Ensamble global
    n_nodes = nodes.shape[0]
    K = assemble_global_matrix(nodes, elements)
    f = np.zeros(n_nodes)

    # 4) Conjuntos de Dirichlet (frontera y placas internas)
    boundary_nodes = []
    boundary_values = []
    for i in range(n_nodes):
        xi, yi = nodes[i]
        if (np.isclose(xi, 0.0) or np.isclose(xi, Lx) or
            np.isclose(yi, 0.0) or np.isclose(yi, Ly)):
            boundary_nodes.append(i); boundary_values.append(0.0)
        elif 2.0 <= xi <= 6.0 and np.isclose(yi, 2.0):
            boundary_nodes.append(i); boundary_values.append(Vp1)
        elif 2.0 <= xi <= 6.0 and np.isclose(yi, 4.0):
            boundary_nodes.append(i); boundary_values.append(Vp2)

    # 5) Resolver con reducción de Dirichlet
    V = solve_with_dirichlet_reduction(K, f, boundary_nodes, boundary_values)

    # 6) Postproceso y guardado como imagen
    V_grid = V.reshape((ny, nx))
    trazar_potencial_y_campo(
        X, Y, V_grid, Vp1, Vp2,
        ruta_salida=r"D:\stuff\Semestre 6\Metodos\Metodos-Numericos-2025-1-Samuel-Gil-Juan-Alvarez-\MIniproyecto2\Punto 1\docs\laplace2d_solucion.png",
        dibujar_placas=False
    )

    # Si quieres exportar datos numéricos:
    # np.savetxt('laplace2d_solution.csv',
    #            np.column_stack((nodes, V)), delimiter=',', header='x,y,V', comments='')


if __name__ == "__main__":
    main()
