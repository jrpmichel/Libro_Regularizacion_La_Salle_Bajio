# Álg. U3 · Sistemas lineales y rango

[![Abrir en Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jrpmichel/Libro_Regularizacion_La_Salle_Bajio/blob/main/alg/u3_sistemas.ipynb)

Notebook completo: [`u3_sistemas.ipynb`](../u3_sistemas.ipynb). Cada celda de código del notebook viene de uno de los archivos `NN_*.py` de esta carpeta, en orden. Los archivos de [`fragmentos/`](fragmentos) son los programas que el libro imprime; cada uno corre solo, sin el notebook.

## Archivos del notebook

| Archivo | Contenido |
|---|---|
| [`00_utilidades.py`](00_utilidades.py) | utilidades compartidas |
| [`01_gauss.py`](01_gauss.py) | sección 3.1 eliminación de Gauss y Gauss-Jordan |
| [`02_cramer.py`](02_cramer.py) | sección 3.2 regla de Cramer y su costo |
| [`03_rango.py`](03_rango.py) | sección 3.3 ecuaciones vectorial y matricial; rango |
| [`04_clasificacion.py`](04_clasificacion.py) | sección 3.4 clasificación de sistemas y sistemas homogéneos |
| [`05_geometria.py`](05_geometria.py) | sección 3.5 interpretación geométrica en 2D y 3D |
| [`99_verificacion_final.py`](99_verificacion_final.py) | verificación final |

## Fragmentos que imprime el libro

| Archivo | Dónde aparece |
|---|---|
| [`alg-u3-l1-prueba.py`](fragmentos/alg-u3-l1-prueba.py) | Laboratorio del Error 3.1 (prueba con codigo) |
| [`alg-u3-l2-prueba.py`](fragmentos/alg-u3-l2-prueba.py) | Laboratorio del Error 3.2 (prueba con codigo) |
| [`alg-u3-l3-prueba.py`](fragmentos/alg-u3-l3-prueba.py) | Laboratorio del Error 3.3 (prueba con codigo) |
| [`alg-u3-r01-verificacion.py`](fragmentos/alg-u3-r01-verificacion.py) | problema resuelto ALG-U3-01 (verificacion) |
| [`alg-u3-r02-verificacion.py`](fragmentos/alg-u3-r02-verificacion.py) | problema resuelto ALG-U3-02 (verificacion) |
| [`alg-u3-r03-verificacion.py`](fragmentos/alg-u3-r03-verificacion.py) | problema resuelto ALG-U3-03 (verificacion) |
| [`alg-u3-r04-verificacion.py`](fragmentos/alg-u3-r04-verificacion.py) | problema resuelto ALG-U3-04 (verificacion) |
| [`alg-u3-r05-verificacion.py`](fragmentos/alg-u3-r05-verificacion.py) | problema resuelto ALG-U3-05 (verificacion) |
| [`alg-u3-s1-gauss.py`](fragmentos/alg-u3-s1-gauss.py) | subtema 3.1: eliminacion de Gauss-Jordan paso a paso |
| [`alg-u3-s2-cramer.py`](fragmentos/alg-u3-s2-cramer.py) | subtema 3.2: regla de Cramer con determinantes |
| [`alg-u3-s3-rango.py`](fragmentos/alg-u3-s3-rango.py) | subtema 3.3: Ax como combinacion de columnas y rango |
| [`alg-u3-s4-clasificacion.py`](fragmentos/alg-u3-s4-clasificacion.py) | subtema 3.4: clasificacion con rangos y solucion general |
| [`alg-u3-s5-planos.py`](fragmentos/alg-u3-s5-planos.py) | subtema 3.5: tres planos y su punto comun |

Para ejecutar un fragmento en tu computadora: `python fragmentos/<archivo>.py`. En Colab, copia su contenido en una celda de código.
