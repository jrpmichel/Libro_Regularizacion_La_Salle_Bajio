# Dif. U5 · Reglas de derivación

[![Abrir en Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jrpmichel/Libro_Regularizacion_La_Salle_Bajio/blob/main/dif/u5_reglas.ipynb)

Notebook completo: [`u5_reglas.ipynb`](../u5_reglas.ipynb). Cada celda de código del notebook viene de uno de los archivos `NN_*.py` de esta carpeta, en orden. Los archivos de [`fragmentos/`](fragmentos) son los programas que el libro imprime; cada uno corre solo, sin el notebook.

## Archivos del notebook

| Archivo | Contenido |
|---|---|
| [`00_utilidades.py`](00_utilidades.py) | utilidades compartidas |
| [`01_potencia.py`](01_potencia.py) | sección 5.1 constante, potencia y suma |
| [`02_producto_cociente.py`](02_producto_cociente.py) | sección 5.2 producto y cociente |
| [`03_cadena.py`](03_cadena.py) | sección 5.3 regla de la cadena |
| [`04_orden_superior.py`](04_orden_superior.py) | sección 5.4 derivadas de orden superior |
| [`05_exp_log_trig.py`](05_exp_log_trig.py) | sección 5.5 exponencial, logaritmo y trigonométricas |
| [`06_implicita.py`](06_implicita.py) | sección 5.6 derivación implícita |
| [`07_logaritmica.py`](07_logaritmica.py) | sección 5.7 derivación logarítmica |
| [`99_verificacion_final.py`](99_verificacion_final.py) | verificación final |

## Fragmentos que imprime el libro

| Archivo | Dónde aparece |
|---|---|
| [`dif-u5-l1-prueba.py`](fragmentos/dif-u5-l1-prueba.py) | Laboratorio del Error 5.1 (prueba con código) |
| [`dif-u5-l2-prueba.py`](fragmentos/dif-u5-l2-prueba.py) | Laboratorio del Error 5.2 (prueba con código) |
| [`dif-u5-l3-prueba.py`](fragmentos/dif-u5-l3-prueba.py) | Laboratorio del Error 5.3 (prueba con código) |
| [`dif-u5-r01-verificacion.py`](fragmentos/dif-u5-r01-verificacion.py) | problema resuelto DIF-U5-01 (verificación) |
| [`dif-u5-r02-verificacion.py`](fragmentos/dif-u5-r02-verificacion.py) | problema resuelto DIF-U5-02 (verificación) |
| [`dif-u5-r03-verificacion.py`](fragmentos/dif-u5-r03-verificacion.py) | problema resuelto DIF-U5-03 (verificación) |
| [`dif-u5-r04-verificacion.py`](fragmentos/dif-u5-r04-verificacion.py) | problema resuelto DIF-U5-04 (verificación) |
| [`dif-u5-r05-verificacion.py`](fragmentos/dif-u5-r05-verificacion.py) | problema resuelto DIF-U5-05 (verificación) |
| [`dif-u5-r06-verificacion.py`](fragmentos/dif-u5-r06-verificacion.py) | problema resuelto DIF-U5-06 (verificación) |
| [`dif-u5-r07-verificacion.py`](fragmentos/dif-u5-r07-verificacion.py) | problema resuelto DIF-U5-07 (verificación) |
| [`dif-u5-s1-potencia.py`](fragmentos/dif-u5-s1-potencia.py) | subtema 5.1: regla de la potencia contra el cociente incremental |
| [`dif-u5-s2-producto.py`](fragmentos/dif-u5-s2-producto.py) | subtema 5.2: regla del producto contra el producto de derivadas |
| [`dif-u5-s3-cadena.py`](fragmentos/dif-u5-s3-cadena.py) | subtema 5.3: la derivada de (x^2 + 1)^10 en x = 1 (semilla) |
| [`dif-u5-s4-orden-superior.py`](fragmentos/dif-u5-s4-orden-superior.py) | subtema 5.4: derivadas sucesivas de un polinomio |
| [`dif-u5-s5-exp-trig.py`](fragmentos/dif-u5-s5-exp-trig.py) | subtema 5.5: pendientes de e^x, 2^x y del seno (radianes y grados) |
| [`dif-u5-s6-implicita.py`](fragmentos/dif-u5-s6-implicita.py) | subtema 5.6: pendiente de la circunferencia x^2 + y^2 = 25 |
| [`dif-u5-s7-logaritmica.py`](fragmentos/dif-u5-s7-logaritmica.py) | subtema 5.7: derivada de x^x por derivación logarítmica |

Para ejecutar un fragmento en tu computadora: `python fragmentos/<archivo>.py`. En Colab, copia su contenido en una celda de código.
