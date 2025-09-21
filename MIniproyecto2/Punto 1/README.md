# Punto 1/ Solución de la ecuación de Laplace
---
Este problema consistió en usar un esquema de diferencias finitas para solucionar la ecuación de Laplace para un capacitor de placas planas paralelas en el vacío.

### Descripción del problema
---
La ecuación en cuestión es:

∇²V = 0

El problema fue modelado en 2 dimensiones, considerando las placas como unidimensionales, de longitud l = 4. Las placas están centras en x = 4, posicionadas en y = 2 y y = 4, y con potenciales Vp1 y Vp2, respectivamente. El dominio del problema fue entonces definido como:

Ω = [0,8] × [0,6]

Se estableció la condición de frontera de Dirichlet:

V(x,y) = 0, ∀(x,y) ∈ ∂Ω

Y se establecieron también las condiciones internas para los potenciales en las placas:

V(x,2) = Vₚ₁, ∀x ∈ [2,6]  
V(x,4) = Vₚ₂, ∀x ∈ [2,6]


### Metodología

---

### Resultados y Análisis

---

### Archivos del Proyecto

* `Laplace2D.py`: Contiene el código fuente que implementa los algoritmos de solución y gráfico de la ecuación.

---

### Instrucciones de Uso

1.  Clona este repositorio.
2.  Instala las dependencias necesarias: `numpy`, `scipy`, `matplotlib`.
3.  Ejecuta el script principal `VolumenChocolatera.py`.

