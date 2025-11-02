# Punto 1/ Solución de la ecuación de Poisson
---

### Descripción del problema
---
El problema consiste en desarrollar un código de Python para resolver la ecuación de Poisson en 1 dimensión, con condiciones de Dirichlet, y usando el Método de Elementos Finitos con elementos lineales.

$$
-\frac{d^2u}{dx^2} = f(x), \quad x \in [0,1], \quad u(0) = u(1) = 0
$$

Esto, con el objetivo de observar el comportamiento de la solución numérica obtenida frente a la solución analítica, empleando una malla no uniforme generada aleatoriamente.

### Metodología
---
El código implementa una formulación del Método de Elementos Finitos (FEM) basada directamente en la forma débil del problema, donde la matriz de rigidez se define como:

$$
K_{ij} = \langle \phi_i', \phi_j' \rangle
$$

y el vector de términos independientes como:

$$
b_i = \langle f, \phi_i \rangle
$$

Para la función base lineal tipo “hat”, la matriz de rigidez global se obtiene a partir de los tamaños de elemento no uniformes $$\(h_i = x_i - x_{i-1}\)$$:

$$
\[
K_{i,i-1} = -\frac{1}{h_i}, \quad
K_{i,i} = \frac{1}{h_i} + \frac{1}{h_{i+1}}, \quad
K_{i,i+1} = -\frac{1}{h_{i+1}}
\]
$$

La malla se genera a partir de una distribución gaussiana centrada en los puntos de una malla uniforme, y los términos del vector \(b\) se calculan mediante **cuadratura gaussiana** en cada elemento.

Finalmente, el sistema lineal \( K u = b \) se resuelve utilizando `numpy.linalg.solve`, y se comparan los resultados con la solución analítica:

$$
u(x) = -\frac{2}{3}x^3 + 3x^2 - \frac{7}{3}x
$$

---

### Resultados y análisis
---
El programa genera dos gráficos principales:

1. **`poisson_1D.png`**: muestra la comparación entre la solución numérica obtenida por FEM y la solución analítica exacta. Se observa una buena concordancia, especialmente al aumentar el número de elementos.
2. **`load_function.png`**: representa la función de carga \( f(x) = 4x - 6 \), utilizada como término fuente en la ecuación de Poisson.


La variación aleatoria en la malla permite observar la robustez del método frente a discretizaciones no uniformes, manteniendo la estabilidad y precisión de la solución.
---

### Archivos del Proyecto

* `script.py`: Contiene el código fuente que implementa los algoritmos de solución y gráficos de las soluciones.
* `poisson_1D.png`: Contiene el gráfico de la solución obtenida con el Método de Elementos Finitos junto con la solución analítica de la ecuación, para hacer una comparación cualitativa.
* `load_function.png`: Contiene el gráfico de la función de carga $$ f(x) $$ usada, graficada en el dominio de la solución.

---

### Instrucciones de Uso

1.  Clona este repositorio.
2.  Instala las dependencias necesarias: `numpy`, `scipy`, `matplotlib`.
3.  Ejecuta el script principal `script.py`.
