import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as st
import pandas as pd
from scipy.spatial import Delaunay
from mpl_toolkits.mplot3d import Axes3D  # Import para gráfico 3D

# ========================================
# Parámetros del problema
# ========================================
R = 0.6 # [m]
alpha = 0.01 # [N/m²]
beta = 50 # [1/m²]
gamma = 0.3 
theta = np.pi / 4  # dirección del máximo cargue (arriba) [rad]

# ========================================
# 1. Definir dominio hexagonal
# ========================================
def define_hexagon(radius):
    points = []
    for i in range(6):
        angle = np.pi / 3 * i
        x = np.sin(angle) * radius
        y = np.cos(angle) * radius
        points.append((x, y))
    return np.array(points)

hex_points = define_hexagon(R)

# ========================================
# 2. Función de carga
# ========================================
def poisson_equation(x, y):
    r = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return alpha * (1 + beta * r**2 + gamma * np.cos(phi - theta))

# ========================================
# 3. Generar puntos internos y malla
# ========================================
def point_inside_polygon(x, y, polygon):
    n = len(polygon)
    inside = False
    for i in range(n):
        x1, y1 = polygon[i]
        x2, y2 = polygon[(i + 1) % n]
        if ((y1 > y) != (y2 > y)) and (x < (x2 - x1) * (y - y1) / (y2 - y1 + 1e-12) + x1):
            inside = not inside
    return inside

def generate_points_in_hexagon(hex_points, n_points=1200):
    min_x, min_y = np.min(hex_points, axis=0)
    max_x, max_y = np.max(hex_points, axis=0)
    pts = []
    while len(pts) < n_points:
        x, y = np.random.uniform(min_x, max_x), np.random.uniform(min_y, max_y)
        if point_inside_polygon(x, y, hex_points):
            pts.append((x, y))
    return np.array(pts)

internal_points = generate_points_in_hexagon(hex_points, n_points=1200)
points = np.vstack([hex_points, internal_points])
triang = Delaunay(points)
vertices = triang.points
elements = triang.simplices

# ========================================
# 4. Ensamblaje FEM
# ========================================
def distance_point_to_segment(p, a, b):
    ap = p - a
    ab = b - a
    t = np.dot(ap, ab) / np.dot(ab, ab)
    t = np.clip(t, 0, 1)
    closest = a + t * ab
    return np.linalg.norm(p - closest)

def fem_solver(vertices, elements, hex_points, tol=0.03):
    
    n_nodes = len(vertices)
    K = np.zeros((n_nodes, n_nodes))
    F = np.zeros(n_nodes)

    for tri in elements:
        coords = vertices[tri]
        x = coords[:, 0]
        y = coords[:, 1]
        area = 0.5 * np.linalg.det([[1, x[0], y[0]],
                                    [1, x[1], y[1]],
                                    [1, x[2], y[2]]])
        if area <= 1e-5:  # evita áreas negativas o degeneradas
            continue
        area = abs(area)

        b = np.array([y[1]-y[2], y[2]-y[0], y[0]-y[1]])
        c = np.array([x[2]-x[1], x[0]-x[2], x[1]-x[0]])
        Ke = (1 / (4 * area)) * (np.outer(b, b) + np.outer(c, c))

        f_val = poisson_equation(np.mean(x), np.mean(y))
        Fe = np.full(3, f_val * area / 3)

        for i in range(3):
            F[tri[i]] += Fe[i]
            for j in range(3):
                K[tri[i], tri[j]] += Ke[i, j]


    # Condición de Dirichlet en todo el borde
    for i, p in enumerate(vertices):
        for j in range(len(hex_points)):
            a = hex_points[j]
            b = hex_points[(j + 1) % len(hex_points)]
            if distance_point_to_segment(p, a, b) < tol:
                K[i, :] = 0
                K[:, i] = 0
                K[i, i] = 1
                F[i] = 0
                break

    U = np.linalg.solve(K, F)
    return U

U = fem_solver(vertices, elements, hex_points)

# ========================================
# 5. Visualización 2D
# ========================================

# --- CARGA ---
fig, ax = plt.subplots(figsize=(7, 6))
X, Y = np.meshgrid(np.linspace(-R, R, 250), np.linspace(-R, R, 250))
Z = poisson_equation(X, Y)
c = ax.pcolormesh(X, Y, Z, cmap='viridis', shading='auto')
ax.plot(*hex_points.T, 'w-', lw=1.5)
ax.set_aspect('equal')
ax.set_title("Distribución de carga $f(x,y)$")
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]")
fig.colorbar(c, ax=ax, label='Carga (N/m²)')
plt.tight_layout()
plt.show()

# --- SOLUCIÓN FEM 2D ---
fig, ax = plt.subplots(figsize=(7, 6))
t = ax.tripcolor(vertices[:, 0], vertices[:, 1], elements, U,
                 shading='gouraud', cmap='inferno')
ax.plot(*hex_points.T, 'w-', lw=1.5)
ax.set_aspect('equal')
ax.set_title('Solución FEM de la ecuación de Poisson\n(u=0 en toda la frontera)')
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
fig.colorbar(t, ax=ax, label='Deformación $u(x,y)$ [m]')
plt.tight_layout()
plt.show()

# ========================================
# 6. Visualización 3D: deformación en metros
# ========================================

fig = plt.figure(figsize=(9, 7))
ax = fig.add_subplot(111, projection='3d')

# Crear la superficie FEM 3D
trisurf = ax.plot_trisurf(vertices[:, 0], vertices[:, 1], U,
                          triangles=elements, cmap='inferno', linewidth=0.1, antialiased=True)

# Escala y rotación
ax.set_title('Superficie 3D de deformación del espejo\nEcuación de Poisson - FEM', pad=15)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('Deformación $u(x,y)$ [m]')
ax.view_init(elev=30, azim=60)

# Colorbar y texto informativo
cbar = fig.colorbar(trisurf, ax=ax, shrink=0.6, pad=0.1)
cbar.set_label('Deformación [m]', rotation=270, labelpad=15)

u_max = np.max(U)
fig.text(0.1, 0.92, f"Deformación máxima: {u_max:.3e} m", color='black', fontsize=11)

plt.tight_layout()
plt.show()

# ========================================
# 7. Análisis Monte Carlo (malla fija, R variable, sin normalización)
# ========================================

def boundary_nodes_mask(vertices, hex_points, tol=0.03):
    mask = np.zeros(len(vertices), dtype=bool)
    for i, p in enumerate(vertices):
        for j in range(len(hex_points)):
            a = hex_points[j]
            b = hex_points[(j + 1) % len(hex_points)]
            if distance_point_to_segment(p, a, b) < tol:
                mask[i] = True
                break
    return mask

bmask = boundary_nodes_mask(vertices, hex_points, tol=0.03)

# --- Matriz K fija ---
def assemble_stiffness(vertices, elements):
    n = len(vertices)
    K = np.zeros((n, n))
    for tri in elements:
        coords = vertices[tri]
        x, y = coords[:, 0], coords[:, 1]
        A = 0.5 * np.linalg.det([[1, x[0], y[0]],
                                 [1, x[1], y[1]],
                                 [1, x[2], y[2]]])
        if A <= 1e-8:
            continue
        A = abs(A)
        b = np.array([y[1]-y[2], y[2]-y[0], y[0]-y[1]])
        c = np.array([x[2]-x[1], x[0]-x[2], x[1]-x[0]])
        Ke = (1 / (4 * A)) * (np.outer(b, b) + np.outer(c, c))
        for i in range(3):
            for j in range(3):
                K[tri[i], tri[j]] += Ke[i, j]
    return K

K = assemble_stiffness(vertices, elements)

Kbc = K.copy()
idx = np.where(bmask)[0]
Kbc[idx, :] = 0
Kbc[:, idx] = 0
Kbc[idx, idx] = 1

# --- Monte Carlo ---
N = 100
u_max_values = []
rng = np.random.default_rng(42)

R_vals     = rng.normal(0.6, 0.02,  N)     # [m]
alpha_vals = rng.normal(0.01, 0.001,    N)     # [N/m²]
beta_vals  = rng.normal(50,  5,     N)     # [1/m²]
gamma_vals = rng.normal(0.3, 0.05,  N)     # adimensional
theta_vals = rng.normal(np.pi/4, np.pi/12, N)    # [rad]

def assemble_load(vertices, elements, R, alpha, beta, gamma, theta):
    F = np.zeros(len(vertices))
    for tri in elements:
        coords = vertices[tri]
        x, y = coords[:, 0], coords[:, 1]
        A = 0.5 * np.linalg.det([[1, x[0], y[0]],
                                 [1, x[1], y[1]],
                                 [1, x[2], y[2]]])
        if A <= 1e-8:
            continue
        A = abs(A)
        xc, yc = np.mean(x), np.mean(y)
        r = np.hypot(xc, yc)
        phi = np.arctan2(yc, xc)
        # SIN normalización
        f_val = alpha * (1 + beta * r**2 + gamma * np.cos(phi - theta))
        Fe = np.full(3, f_val * A / 3)
        for i in range(3):
            F[tri[i]] += Fe[i]
    return F

for i in range(N):
    R_i, a_i, b_i, g_i, th_i = R_vals[i], alpha_vals[i], beta_vals[i], gamma_vals[i], theta_vals[i]
    F = assemble_load(vertices, elements, R_i, a_i, b_i, g_i, th_i)
    F[bmask] = 0.0
    try:
        U_i = np.linalg.solve(Kbc, F)
        u_max_values.append(np.max(U_i))
    except np.linalg.LinAlgError:
        continue
    print(f"Simulación {i+1}/{N} completada")

u_max_values = np.array(u_max_values)

# --- Estadística ---
mean_u = np.mean(u_max_values)
std_u  = np.std(u_max_values, ddof=1)
ci95   = st.t.interval(0.95, len(u_max_values)-1, loc=mean_u, scale=std_u/np.sqrt(len(u_max_values)))

print("\n===== RESULTADOS MONTE CARLO =====")
print(f"Simulaciones válidas: {len(u_max_values)} / {N}")
print(f"Media de u_max = {mean_u:.4e} m")
print(f"Desviación estándar = {std_u:.4e} m")
print(f"IC 95% = ({ci95[0]:.4e}, {ci95[1]:.4e}) m")

# --- Histograma ---
plt.figure(figsize=(8,5))
sns.histplot(u_max_values, bins=15, color='orange', edgecolor='k')
plt.axvline(mean_u, color='r', linestyle='--', label=f"Media = {mean_u:.3e} m")
plt.title("Distribución Monte Carlo de la deformación máxima (u_max)")
plt.xlabel("u_max [m]")
plt.ylabel("Frecuencia")
plt.legend()
plt.tight_layout()
plt.show()

# ========================================
# 8. Análisis de correlación entre parámetros y deformación máxima
# ========================================

# Crear DataFrame con parámetros y resultados
df_corr = pd.DataFrame({
    "R": R_vals[:len(u_max_values)],
    "alpha": alpha_vals[:len(u_max_values)],
    "beta": beta_vals[:len(u_max_values)],
    "gamma": gamma_vals[:len(u_max_values)],
    "theta": theta_vals[:len(u_max_values)],
    "u_max": u_max_values
})

# Calcular matriz de correlación
corr_matrix = df_corr.corr(method="pearson")

# Extraer correlación de u_max con los parámetros
corr_u = corr_matrix["u_max"].drop("u_max").sort_values(ascending=False)

print("\n===== CORRELACIÓN ENTRE PARÁMETROS Y u_max =====")
print(corr_u)

# --- Visualización ---
plt.figure(figsize=(7,5))
sns.barplot(x=corr_u.values, y=corr_u.index, hue=corr_u.index, palette="viridis", legend=False)
plt.title("Correlación entre parámetros y deformación máxima (u_max)")
plt.xlabel("Coeficiente de correlación de Pearson")
plt.ylabel("Parámetro")
plt.grid(axis="x", linestyle="--", alpha=0.6)
plt.tight_layout()
plt.show()

# ========================================
# 9. Estudio Paramétrico para todos los parámetros (α, β, γ, R, θ)
# ========================================

def estudio_parametrico(nombre_param, valores, R_mean=0.6, alpha_mean=100, beta_mean=50, gamma_mean=0.3, theta_mean=np.pi/4):
    """Ejecuta el FEM variando un único parámetro y retorna los resultados."""
    u_max_list = []
    for i, val in enumerate(valores):
        F = np.zeros(len(vertices))
        for tri in elements:
            coords = vertices[tri]
            x, y = coords[:, 0], coords[:, 1]
            A = 0.5 * np.linalg.det([[1, x[0], y[0]],
                                     [1, x[1], y[1]],
                                     [1, x[2], y[2]]])
            if A <= 1e-8:
                continue
            A = abs(A)
            xc, yc = np.mean(x), np.mean(y)
            r = np.hypot(xc, yc)
            phi = np.arctan2(yc, xc)

            # Selección de parámetro variable
            if nombre_param == "alpha":
                f_val = val * (1 + beta_mean * r**2 + gamma_mean * np.cos(phi - theta_mean))
            elif nombre_param == "beta":
                f_val = alpha_mean * (1 + val * r**2 + gamma_mean * np.cos(phi - theta_mean))
            elif nombre_param == "gamma":
                f_val = alpha_mean * (1 + beta_mean * r**2 + val * np.cos(phi - theta_mean))
            elif nombre_param == "R":
                f_val = alpha_mean * (1 + beta_mean * (r**2 / val) + gamma_mean * np.cos(phi - theta_mean))
            elif nombre_param == "theta":
                f_val = alpha_mean * (1 + beta_mean * r**2 + gamma_mean * np.cos(phi - val))
            else:
                raise ValueError("Parámetro no reconocido")

            Fe = np.full(3, f_val * A / 3)
            for k in range(3):
                F[tri[k]] += Fe[k]
        F[bmask] = 0.0

        try:
            U_i = np.linalg.solve(Kbc, F)
            u_max_list.append(np.max(U_i))
        except np.linalg.LinAlgError:
            u_max_list.append(np.nan)

        print(f"Simulación {nombre_param} {i+1}/{len(valores)} completada")
    return np.array(u_max_list)

# --- Ejecutar estudios ---
u_alpha = estudio_parametrico("alpha", alpha_vals)
u_beta  = estudio_parametrico("beta",  beta_vals)
u_gamma = estudio_parametrico("gamma", gamma_vals)
u_R     = estudio_parametrico("R",     R_vals)
u_theta = estudio_parametrico("theta", theta_vals)

# --- Graficar resultados comparativos ---
parametros = {
    "alpha": (alpha_vals, u_alpha, "α [N/m²]", True),
    "beta":  (beta_vals,  u_beta,  "β [1/m²]", True),
    "gamma": (gamma_vals, u_gamma, "γ [-]", True),
    "R":     (R_vals,     u_R,     "R [m]", True),
    "theta": (theta_vals, u_theta, "θ [rad]", False)  # sin regresión lineal
}

fig, axs = plt.subplots(3, 2, figsize=(11, 12))
axs = axs.flatten()

for i, (param, (x, y, label, lineal)) in enumerate(parametros.items()):
    sns.scatterplot(x=x, y=y, ax=axs[i], color="teal", s=40, alpha=0.7, edgecolor="k", linewidth=0.3)
    if lineal:
        sns.regplot(x=x, y=y, ax=axs[i], scatter=False, color="red", ci=None)
    axs[i].set_title(f"{param} vs $u_{{max}}$")
    axs[i].set_xlabel(label)
    axs[i].set_ylabel("$u_{max}$ [m]")
    axs[i].grid(True, linestyle="--", alpha=0.5)

# Eliminar subplot vacío
for j in range(len(parametros), len(axs)):
    fig.delaxes(axs[j])

plt.suptitle("Estudio Paramétrico Univariado – Influencia de cada parámetro en la deformación máxima", fontsize=13)
plt.tight_layout(rect=[0, 0, 1, 0.97])
plt.show()
