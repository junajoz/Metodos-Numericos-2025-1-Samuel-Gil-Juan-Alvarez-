# Punto 1/ Solución de la ecuación de Laplace
---
Este problema consistió en usar un esquema de diferencias finitas para solucionar la ecuación de Laplace para un capacitor de placas planas paralelas en el vacío.

### Descripción del problema
---
La ecuación en cuestión es:

∇²V = 0

El problema fue modelado en 2 dimensiones, considerando las placas como unidimensionales, de longitud l = 4. Las placas están centras en x = 4, posicionadas en y = 2 y y = 4, y con potenciales Vₚ₁ y Vₚ₂, respectivamente. El dominio del problema fue entonces definido como:

Ω = [0,8] × [0,6]

Se estableció la condición de frontera de Dirichlet:

V(x,y) = 0, ∀(x,y) ∈ ∂Ω

Y se establecieron también las condiciones internas para los potenciales en las placas:

V(x,2) = Vₚ₁, ∀x ∈ [2,6]  
V(x,4) = Vₚ₂, ∀x ∈ [2,6]


### Metodología
---
Para resolver la ecuación de Laplace en dos dimensiones se empleó un esquema de aproximación numérica basado en el método de diferencias finitas, implementado mediante un mallado rectangular del dominio y la posterior discretización en elementos triangulares lineales (P1). Esto permite transformar la ecuación diferencial en un sistema lineal de ecuaciones que puede resolverse de manera matricial.

El código se organiza de la siguiente forma:

1. **Generación de la malla**: el dominio Ω = [0,8] × [0,6] se discretiza en una malla rectangular de puntos, y cada celda rectangular se divide en dos triángulos para la formulación numérica.  
2. **Ensamble de la matriz global**: se calcula la matriz de rigidez local para cada triángulo y se ensambla en una matriz global dispersa.  
3. **Condiciones de frontera**: se imponen condiciones de Dirichlet, asignando V = 0 en el contorno del dominio y valores fijos Vₚ₁ y Vₚ₂ en las placas internas situadas en y = 2 y y = 4, respectivamente.  
4. **Resolución del sistema**: el sistema lineal reducido se resuelve para obtener el potencial en cada nodo libre de la malla.  
5. **Postprocesado y visualización**: se reconstruye el campo potencial sobre la malla y se grafican tanto las líneas equipotenciales como el campo eléctrico asociado, representando el comportamiento físico de un capacitor de placas paralelas.

Este enfoque permite aproximar la solución de la ecuación de Laplace en el dominio definido y visualizar tanto la distribución del potencial eléctrico como el campo resultante.

---

### Resultados y Análisis

---

### Archivos del Proyecto

* `Laplace2D.py`: Contiene el código fuente que implementa los algoritmos de solución y gráfico de la ecuación.
* `laplace2d_solucion.png`: Contiene la imagen representando la solución del problema como un mapa de calor, donde el color indica la magnitud del potencial eléctrico en un punto del plano xy.

---

### Instrucciones de Uso

1.  Clona este repositorio.
2.  Instala las dependencias necesarias: `numpy`, `scipy`, `matplotlib`.
3.  Ejecuta el script principal `VolumenChocolatera.py`.

