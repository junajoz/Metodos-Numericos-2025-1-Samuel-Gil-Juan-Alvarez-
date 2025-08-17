import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simpson

# -------------------------------
# Datos originales
# -------------------------------
x = np.array([0., 3.32994924, 8.53299492, 17.4822335, 24.76649746,
              30.3857868, 32.46700507, 36.42131979, 48.49238578, 57.44162436,
              66.18274111, 66.8071066, 73.88324873, 76.58883248, 86.78680203,
              98.23350253, 106.76649746, 114.6751269, 123.62436548, 134.4467005,
              144.02030456, 151.51269035, 158.58883248, 162.33502537, 165.66497461,
              169.61928933, 173.1573604, 176.27918781, 178.15228425, 180.0253807,
              181.69035532, 182.73096446, 183.35532994, 183.56345177])

y = np.array([61.39593908, 59.52284264, 57.02538071, 53.9035533,  52.23857868, 51.40609137,
              51.19796954, 50.98984771, 51.19796954, 52.65482233, 55.15228426, 55.36040609,
              57.85786802, 58.89847715, 62.85279187, 65.3502538,  66.18274111, 66.39086294,
              65.97461929, 64.30964467, 61.39593908, 57.85786802, 53.48730964, 50.57360406,
              47.65989847, 43.70558375, 39.33502538, 34.34010152, 30.59390863, 25.59898477,
              19.35532995, 13.52791878,  7.70050761,  0.])

# -------------------------------
# Ordenar y limpiar datos
# -------------------------------
data = np.column_stack((x, y))
data = data[np.argsort(data[:, 0])]
_, idx = np.unique(data[:, 0], return_index=True)
data = data[idx]
x_sorted = data[:, 0]
y_sorted = data[:, 1]

# -------------------------------
# Interpolación por tramos grado 2
# -------------------------------
x_interp = []
y_interp = []
piecewise_polys = []

for i in range(len(x_sorted) - 2):
    xi = x_sorted[i:i+3]
    yi = y_sorted[i:i+3]
    coeffs = np.polyfit(xi, yi, 2)
    poly = np.poly1d(coeffs)
    piecewise_polys.append((poly, xi[0], xi[-1]))

    xi_dense = np.linspace(xi[0], xi[-1], 50)
    yi_dense = poly(xi_dense)

    x_interp.extend(xi_dense)
    y_interp.extend(yi_dense)

x_interp = np.array(x_interp)
y_interp = np.array(y_interp)

# -------------------------------
# Calcular volumen de revolución
# -------------------------------
V_mm3 = np.pi * simpson(y_interp**2, x_interp)  # mm³
V_mL = V_mm3 / 1000  # mL

# -------------------------------
# Resultados
# -------------------------------

print("Polinomios por tramo (grado 2):")
for i, (poly, x_start, x_end) in enumerate(piecewise_polys):
    print(f"Tramo {i+1} ({x_start:.3f} a {x_end:.3f}): y = {poly}")
error = abs((V_mL-2200)/2200)*100

# -------------------------------
# Gráfica
# -------------------------------
plt.figure(figsize=(12, 6))
plt.plot(x_sorted, y_sorted, 'o', label='Datos originales ordenados', markersize=6)
plt.plot(x_interp, y_interp, '-', label='Interpolación por tramos (grado 2)')
plt.xlabel('x (mm)')
plt.ylabel('y (mm)')
plt.title('Interpolación de Lagrange por tramos (grado 2) - Datos limpios')
plt.legend()
plt.grid(True)
plt.show()

print(f"Volumen por revolución ≈ {V_mL:.2f} mL\n")
print(f"Se obtuvo un error de = {error:.2f}%\n")
