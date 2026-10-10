# Álg. U2 · Matrices y determinantes

[![Abrir en Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jrpmichel/Libro_Regularizacion_La_Salle_Bajio/blob/main/alg/u2_matrices.ipynb)

Notebook completo: [`u2_matrices.ipynb`](../u2_matrices.ipynb). Cada celda de código del notebook viene de uno de los archivos `NN_*.py` de esta carpeta, en orden. Los archivos de [`fragmentos/`](fragmentos) son los programas que el libro imprime; cada uno corre solo, sin el notebook.

## Archivos del notebook

| Archivo | Contenido |
|---|---|
| [`00_utilidades.py`](00_utilidades.py) | utilidades compartidas |
| [`01_operaciones.py`](01_operaciones.py) | sección 2.1 operaciones con matrices y transpuesta |
| [`02_determinantes.py`](02_determinantes.py) | sección 2.2 determinantes por Sarrus y por cofactores |
| [`03_inversa.py`](03_inversa.py) | sección 2.3 propiedades del determinante y matriz inversa |
| [`04_geometria.py`](04_geometria.py) | sección 2.4 interpretación geométrica del determinante |
| [`05_lu.py`](05_lu.py) | sección 2.5 factorización LU |
| [`99_verificacion_final.py`](99_verificacion_final.py) | verificación final |

## Fragmentos que imprime el libro

| Archivo | Dónde aparece |
|---|---|
| [`alg-u2-l1-prueba.py`](fragmentos/alg-u2-l1-prueba.py) | Laboratorio del Error 2.1 (prueba con codigo) |
| [`alg-u2-l2-prueba.py`](fragmentos/alg-u2-l2-prueba.py) | Laboratorio del Error 2.2 (prueba con codigo) |
| [`alg-u2-l3-prueba.py`](fragmentos/alg-u2-l3-prueba.py) | Laboratorio del Error 2.3 (prueba con codigo) |
| [`alg-u2-r01-verificacion.py`](fragmentos/alg-u2-r01-verificacion.py) | problema resuelto ALG-U2-01 (verificacion) |
| [`alg-u2-r02-verificacion.py`](fragmentos/alg-u2-r02-verificacion.py) | problema resuelto ALG-U2-02 (verificacion) |
| [`alg-u2-r03-verificacion.py`](fragmentos/alg-u2-r03-verificacion.py) | problema resuelto ALG-U2-03 (verificacion) |
| [`alg-u2-r04-verificacion.py`](fragmentos/alg-u2-r04-verificacion.py) | problema resuelto ALG-U2-04 (verificacion) |
| [`alg-u2-r05-verificacion.py`](fragmentos/alg-u2-r05-verificacion.py) | problema resuelto ALG-U2-05 (verificacion) |
| [`alg-u2-s1-producto.py`](fragmentos/alg-u2-s1-producto.py) | subtema 2.1: producto de matrices, orden de los factores y transpuesta |
| [`alg-u2-s2-determinante.py`](fragmentos/alg-u2-s2-determinante.py) | subtema 2.2: determinante por Sarrus y por cofactores |
| [`alg-u2-s3-inversa.py`](fragmentos/alg-u2-s3-inversa.py) | subtema 2.3: adjunta, inversa y comprobacion A A^-1 = I |
| [`alg-u2-s4-geometria.py`](fragmentos/alg-u2-s4-geometria.py) | subtema 2.4: imagen del cuadrado unitario y area = |det A| |
| [`alg-u2-s5-lu.py`](fragmentos/alg-u2-s5-lu.py) | subtema 2.5: factorizacion LU sin intercambio de renglones |

Para ejecutar un fragmento en tu computadora: `python fragmentos/<archivo>.py`. En Colab, copia su contenido en una celda de código.
