import numpy as np
import matplotlib.pyplot as plt

# Parámetros
a = 2.0
I = 1.0
T = 5.0
dt_values = [1.25, 0.625, 0.3125, 0.15625, 0.078125, 0.0390625]  # Varios Δt

# Solución analítica
t_analytic = np.linspace(0, T, 200)
u_analytic = I * np.exp(-a * t_analytic)

# Solución numérica
methods = ['Forward Euler', 'Crank-Nicolson', 'Backward Euler']
solutions = [[0]*6, [0]*6, [0]*6]

m = 0
for method in methods:  # Ciclo sobre los diferentes métodos
    fig, ax = plt.subplots(1, 2, figsize=(8, 5))
    ax[0].plot(t_analytic, u_analytic, 'k--', label='Solución analítica')

    E_values = np.zeros(len(dt_values))
    for i in range(len(dt_values)):  # Ciclo sobre los diferentes valores de Δt
        dt = dt_values[i]
        N = int(T/dt)
        t = np.linspace(0, N*dt, N + 1)
        u = np.zeros(N + 1)
        u[0] = I

        for n in range(N):  # Ciclo sobre cada punto de la solución
            if method == 'Forward Euler':
                u[n+1] = u[n] - a * dt * u[n]
            elif method == 'Crank-Nicolson':
                u[n+1] = ((1 - 0.5*a*dt) / (1 + 0.5*a*dt)) * u[n]
            elif method == 'Backward Euler':
                u[n+1] = (1 / (1 + a*dt)) * u[n]

        # Guardar soluciones para más adelante
        solutions[m][i] = u
        if i <= 3:
            ax[0].plot(t, u, marker='o', label=f'Δt = {dt}')

        # Análisis del error de cada método
        e = np.zeros(N + 1)
        cumulative_e = np.zeros(N + 2)
        for n in range(N + 1):
            e[n] = I*np.exp(-a*n*dt) - solutions[m][i][n]
            cumulative_e[n + 1] = cumulative_e[n] + e[n]**2
        E_values[i] = np.sqrt(dt*cumulative_e[-1])

    # Encontrar la pendiente de la parte lineal
    x = np.log(dt_values)
    y = np.log(E_values)
    slope, b = np.polyfit(x[1:], y[1:], 1)
    plt.plot(x, slope*x + b)
    ax[1].scatter(x, y)
    print("slope = ", slope)

    ax[1].set_title("Análisis del error")
    ax[1].set_xlabel("Ln(Δt)")
    ax[1].set_ylabel("Ln(E)")
    ax[1].grid(True)
    ax[0].set_xlabel('Tiempo t')
    ax[0].set_ylabel('u(t)')
    ax[0].set_title(f'{method} - Estabilidad')
    ax[0].grid(True)
    ax[0].legend()

    # Anotaciones
    axis = fig.gca()  # Get current axis

    texto = (
        r"$u'(t)=-a\,u(t)$" "\n"
        r"$\frac{2}{a}=1.0$" "\n"
        r"$\frac{1}{a}=0.5$"
    )

    axis.text(
        0.4, 0.95, texto,  # esquina superior derecha del área del eje
        transform=axis.transAxes,  # relativo al eje (0-1)
        ha="right", va="top",
        fontsize=12,
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor="none"),
    )
    fig.tight_layout()
    plt.show()
    m += 1
