# Punto 3: Razón entre Coordenadas y Jacobiano

Este proyecto se centra en la explicación del cambio de variables en integrales múltiples, específicamente la razón entre las coordenadas, a través del uso del **Jacobiano**. El objetivo es proporcionar una explicación contundente, no formal, sobre la relación:

$$dx dy = |J| d\xi d\eta$$

---

### Instrucciones del Proyecto

Este proyecto se basa en un enfoque de lenguaje natural y requiere que se sigan las siguientes instrucciones:

1.  Crea un archivo `.tex` que solo contenga el *prompt* que le darás a un modelo de lenguaje natural para pedirle la tarea. Realiza un *commit* con este *prompt* en el repositorio.
2.  Genera el archivo `.tex` con el modelo de lenguaje y pégalo tal como se lo entrega en el archivo `.tex` anterior, y haz el *commit* ahí.
3.  Realiza todas las ediciones pertinentes, entre las cuales se puede incluir peticiones adicionales al modelo, ediciones manuales, entre otras. Realiza un tercer *commit* con la versión final del `.tex` de tal forma que se observen los cambios entre la versión original y la editada. Además, agrega el PDF generado y menciona el modelo de lenguaje usado.

---

### Prompt Utilizado (Chat GPT)


```text
[PEGA AQUÍ EL PROMPT QUE USaste]
```

### Resumen

La explicación general se enfoca en el concepto del Jacobiano como un factor de escala. Cuando transformas un sistema de coordenadas (como de (x, y) a (ξ, η)) para simplificar una integral, estás cambiando el "tamaño" de los elementos infinitesimales de área (dxdy). El Jacobiano es la matriz de las derivadas parciales que te permite calcular exactamente por cuánto se "estira" o "comprime" el área en esa transformación. Su determinante, ∣J∣, actúa como ese factor de corrección.

### Archivos

En la carpeta /docs se encuentra el .tex donde se desarrolla el punto.
