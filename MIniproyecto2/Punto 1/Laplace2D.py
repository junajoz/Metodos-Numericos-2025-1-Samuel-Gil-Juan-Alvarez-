# -*- coding: utf-8 -*-
"""
Laplace 2D (∇²V=0) en Ω=[0,8]x[0,6] con diferencias finitas centrales (5 puntos).
Condiciones de Dirichlet:
  - V=0 en la frontera ∂Ω
  - Placa inferior:  y=2,  x∈[2,6] → Vp1
  - Placa superior:  y=4,  x∈[2,6] → Vp2
Autores: Juan Álvarez y Samuel Gil
"""

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib.pyplot as plt


def mesh_rectangular(Lx, Ly, nx, ny):
    """Malla cartesiana regular y vectores de coordenadas."""
    x = np.linspace(0.0, Lx, nx)
    y = np.linspace(0.0, Ly, ny)
    X, Y = np.meshgrid(x, y)  # filas→y, cols→x
    return x, y, X, Y, nx, ny

def assemble_A_b(nx, ny, dx, dy, x, y, Vp1, Vp2,
                 x_min_plate=2.0, x_max_plate=6.0,
                 y_plate_vals=(2.0, 4.0)):
    N = nx * ny
    A = sp.lil_matrix((N, N), dtype=float)
    b = np.zeros(N, dtype=float)

    inv_dx2 = 1.0 / (dx * dx)                 # 1/dx^2
    inv_dy2 = 1.0 / (dy * dy)                 # 1/dy^2
    c0 = -2.0 * (inv_dx2 + inv_dy2)           # -2(1/dx^2 + 1/dy^2)

    def pidx(i, j): return j * nx + i

    # Ubicar filas j de las placas (y=2 y y=4)
    placas_j = []
    for yval in y_plate_vals:
        hits = np.where(np.isclose(y, yval))[0]
        placas_j.append(int(hits[0]) if hits.size == 1 else None)

    j_y2, j_y4 = placas_j

    # Recorremos todos los nodos
    for j in range(ny):                       
        for i in range(nx):                   
            p = pidx(i, j)                    

            # --- Frontera ---
            on_left   = (i == 0)              
            on_right  = (i == nx - 1)         
            on_bottom = (j == 0)              
            on_top    = (j == ny - 1)         

            if on_left or on_right or on_bottom or on_top:
                A[p, p] = 1.0                 # Dirichlet frontera
                b[p] = 0.0                    
                continue                      

            # --- Placas internas ---
            es_placa = False                  
            if j_y2 is not None and j == j_y2 and (x_min_plate - 1e-12) <= x[i] <= (x_max_plate + 1e-12):
                A[p, p] = 1.0                 # ← Dirichlet placa y=2
                b[p] = Vp1                    
                es_placa = True               
            if j_y4 is not None and j == j_y4 and (x_min_plate - 1e-12) <= x[i] <= (x_max_plate + 1e-12):
                A[p, p] = 1.0                 # ← Dirichlet placa y=4
                b[p] = Vp2                    
                es_placa = True               

            if es_placa:                      
                continue                      

            # --- Nodo interior: stencil 5 puntos ---
            A[p, p]       = c0                # ← diagonal: -2(1/dx^2+1/dy^2)
            A[p, p - 1]   = inv_dx2           # ← vecino i-1, j
            A[p, p + 1]   = inv_dx2           # ← vecino i+1, j
            A[p, p - nx]  = inv_dy2           # ← vecino i, j-1
            A[p, p + nx]  = inv_dy2           # ← vecino i, j+1
    return A.tocsr(), b

def trazar_potencial_y_campo(X, Y, V_grid, Vp1, Vp2,
                             ruta_salida="laplace2d_solucion.png",
                             dibujar_placas=True):
    """Mapa de potencial + equipotenciales + líneas de campo, y guarda PNG."""
    dx = X[0, 1] - X[0, 0]
    dy = Y[1, 0] - Y[0, 0]
    dV_dy, dV_dx = np.gradient(V_grid, dy, dx)
    Ex = -dV_dx
    Ey = -dV_dy

    fig, ax = plt.subplots(figsize=(8, 6))
    cf = ax.contourf(X, Y, V_grid, levels=50, cmap="jet")
    ax.contour(X, Y, V_grid, levels=18, linewidths=0.3, colors="k", alpha=0.5)

    sp_ = ax.streamplot(X, Y, Ex, Ey, density=1.2, linewidth=0.6,
                        color="k", arrowsize=0.6)
    sp_.lines.set_alpha(0.5)
    sp_.arrows.set_alpha(0.5)

    if dibujar_placas:
        ax.plot([2, 6], [2, 2], linewidth=3, color="w")   # Vp1
        ax.plot([2, 6], [4, 4], linewidth=3, color="w")   # Vp2

    cbar = fig.colorbar(cf, ax=ax)
    cbar.set_label("Potencial [V]")

    ax.set_title("Laplace 2D (FDM): Potencial y líneas de campo")
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_aspect("equal")
    ax.set_xlim(X.min(), X.max())
    ax.set_ylim(Y.min(), Y.max())
    plt.tight_layout()
    plt.savefig(ruta_salida, dpi=300, bbox_inches="tight")
    plt.show()
    plt.close(fig)


def main():
    # Dominio y malla 
    Lx, Ly = 8.0, 6.0
    nx, ny = 81, 61
    Vp1, Vp2 = 1.0, -1.0

    x, y, X, Y, nx, ny = mesh_rectangular(Lx, Ly, nx, ny)
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    # Ensamble directo A y b (incluye bordes y placas)
    A, b = assemble_A_b(nx, ny, dx, dy, x, y, Vp1, Vp2)

    # Resolver A V = b
    V = spla.spsolve(A, b)

    # Post-proceso
    V_grid = V.reshape((ny, nx))
    trazar_potencial_y_campo(
        X, Y, V_grid, Vp1, Vp2,
        ruta_salida=r"D:\stuff\Semestre 6\Metodos\Metodos-Numericos-2025-1-Samuel-Gil-Juan-Alvarez-\MIniproyecto2\Punto 1\docs\laplace2d_solucion.png",
        dibujar_placas=True
    )


if __name__ == "__main__":
    main()
