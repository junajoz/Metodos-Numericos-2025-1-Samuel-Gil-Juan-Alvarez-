# Punto 1/ Cálculo de Volumen por Interpolación Numérica

Este proyecto se enfoca en el cálculo del volumen de un objeto tridimensional con una forma compleja (una chocolatera) utilizando la combinación de dos métodos numéricos: **interpolación de Lagrange** y la **integración por el método trapezoidal/regla de simpson**. El objetivo fue modelar la geometría del objeto a partir de una digitalización 2D para, posteriormente, calcular su volumen.

---

### Metodología

El proceso se ejecutó en los siguientes pasos secuenciales:

1.  **Digitalización de la Geometría**: A partir de una fotografía, se extrajo un conjunto de puntos de datos de manera manual que representaban el perfil de la chocolatera. Estos puntos (`(x, y)`) sirvieron como la entrada de datos para el análisis.

2.  **Interpolación Polinómica**: Para obtener una función continua que describiera la forma del objeto, se aplicó el **método de interpolación de Lagrange**. Este método generó un polinomio único de grado $n-1$ que pasa exactamente por los $n$ puntos de datos tomando subintervalos, para este caso, elegimos parabolas y subdividimos los intervalos en 3 puntos, creando así una curva a trozos que modela el perfil de la chocolatera.

3.  **Cálculo de Volumen por Integración**: El volumen del sólido de revolución se calculó integrando la función obtenida en el paso anterior. Para ello, se empleó el **método del trapecio** y el **método parabolobico**, dos técnicas de integración numérica que aproximan el área bajo la curva dividiéndola en una serie de trapezoides y parabolas para acabar sumando sus áreas. 

    La fórmula utilizada para el volumen de un sólido de revolución es:
    $$V = \int_{a}^{b} \pi [f(x)]^2 dx$$

---

### Resultados y Análisis

El volumen calculado numéricamente se comparó con el volumen real del objeto, resultando en un **error porcentual del 15.19%**.

Este error, considerablemente alto, se atribuye principalmente a las limitaciones en la etapa de digitalización. La naturaleza manual de la extracción de puntos a partir de una imagen bidimensional introduce imprecisiones inherentes que son susceptibles a la dispersión. Pequeñas variaciones en la perspectiva de la fotografía o en la selección de los puntos pueden alterar significativamente la forma del polinomio de Lagrange o su ubicación con respecto al eje coordenado, afectando la precisión del volumen final. Este resultado subraya la sensibilidad de los métodos numéricos a la calidad de los datos de entrada.

---

### Archivos del Proyecto

* `VolumenChocolatera.py`: Contiene el código fuente que implementa los algoritmos de interpolación y cálculo de volumen.
* `ImagenChocolatera.py`: Almacena el código con el que se extrajeron los puntos de la chocolatera.
* `documentos/`: Aquí se encuentra la imagen de la grafica de la curva interpolada y la chocolatera.

---

### Instrucciones de Uso

1.  Clona este repositorio.
2.  Instala las dependencias necesarias: `numpy`, `scipy`, `matplotlib`.
3.  Ejecuta el script principal `VolumenChocolatera.py`.
