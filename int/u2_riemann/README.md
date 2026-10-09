# Int. U2 · Sumas de Riemann

[![Abrir en Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jrpmichel/Libro_Regularizacion_La_Salle_Bajio/blob/main/int/u2_riemann.ipynb)

Notebook completo: [`u2_riemann.ipynb`](../u2_riemann.ipynb). Cada celda de código del notebook viene de uno de los archivos `NN_*.py` de esta carpeta, en orden. Los archivos de [`fragmentos/`](fragmentos) son los programas que el libro imprime; cada uno corre solo, sin el notebook.

## Archivos del notebook

| Archivo | Contenido |
|---|---|
| [`00_utilidades.py`](00_utilidades.py) | utilidades compartidas |
| [`01_sumas.py`](01_sumas.py) | sección 2.1 sumas izquierda, derecha y de punto medio |
| [`02_sigma.py`](02_sigma.py) | sección 2.2 notación sigma |
| [`03_definicion.py`](03_definicion.py) | sección 2.3 la integral definida como límite de sumas |
| [`04_propiedades.py`](04_propiedades.py) | sección 2.4 propiedades de la integral definida |
| [`05_trapecio_simpson.py`](05_trapecio_simpson.py) | sección 2.5 trapecio, Simpson y tabla de error |
| [`99_verificacion_final.py`](99_verificacion_final.py) | verificación final |

## Fragmentos que imprime el libro

| Archivo | Dónde aparece |
|---|---|
| [`int-u2-l1-prueba.py`](fragmentos/int-u2-l1-prueba.py) | Laboratorio del Error 2.1 (prueba con código) |
| [`int-u2-l2-prueba.py`](fragmentos/int-u2-l2-prueba.py) | Laboratorio del Error 2.2 (prueba con código) |
| [`int-u2-l3-prueba.py`](fragmentos/int-u2-l3-prueba.py) | Laboratorio del Error 2.3 (prueba con código) |
| [`int-u2-r01-verificacion.py`](fragmentos/int-u2-r01-verificacion.py) | problema resuelto INT-U2-01 (verificación) |
| [`int-u2-r02-verificacion.py`](fragmentos/int-u2-r02-verificacion.py) | problema resuelto INT-U2-02 (verificación) |
| [`int-u2-r03-verificacion.py`](fragmentos/int-u2-r03-verificacion.py) | problema resuelto INT-U2-03 (verificación) |
| [`int-u2-r04-verificacion.py`](fragmentos/int-u2-r04-verificacion.py) | problema resuelto INT-U2-04 (verificación) |
| [`int-u2-s1-izq-der-med.py`](fragmentos/int-u2-s1-izq-der-med.py) | subtema 2.1: sumas izquierda, derecha y de punto medio |
| [`int-u2-s2-sumas-cerradas.py`](fragmentos/int-u2-s2-sumas-cerradas.py) | subtema 2.2: fórmulas cerradas de las sumas |
| [`int-u2-s3-limite-sumas.py`](fragmentos/int-u2-s3-limite-sumas.py) | subtema 2.3: la integral de x^2 como límite de sumas derechas |
| [`int-u2-s4-propiedades.py`](fragmentos/int-u2-s4-propiedades.py) | subtema 2.4: propiedades de la integral comprobadas con sumas |
| [`int-u2-s5-tabla-error.py`](fragmentos/int-u2-s5-tabla-error.py) | subtema 2.5: tabla de error contra n de cuatro reglas |

Para ejecutar un fragmento en tu computadora: `python fragmentos/<archivo>.py`. En Colab, copia su contenido en una celda de código.
