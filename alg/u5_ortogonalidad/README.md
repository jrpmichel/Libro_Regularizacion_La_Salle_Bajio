# Álg. U5 · Ortogonalidad y mínimos cuadrados

[![Abrir en Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jrpmichel/Libro_Regularizacion_La_Salle_Bajio/blob/main/alg/u5_ortogonalidad.ipynb)

Notebook completo: [`u5_ortogonalidad.ipynb`](../u5_ortogonalidad.ipynb). Cada celda de código del notebook viene de uno de los archivos `NN_*.py` de esta carpeta, en orden. Los archivos de [`fragmentos/`](fragmentos) son los programas que el libro imprime; cada uno corre solo, sin el notebook.

## Archivos del notebook

| Archivo | Contenido |
|---|---|
| [`00_utilidades.py`](00_utilidades.py) | utilidades compartidas |
| [`01_producto_interior.py`](01_producto_interior.py) | sección 5.1 producto interior, longitud y ortogonalidad |
| [`02_proyeccion.py`](02_proyeccion.py) | sección 5.2 proyección ortogonal sobre una recta y sobre un subespacio |
| [`03_gram_schmidt.py`](03_gram_schmidt.py) | sección 5.3 proceso de Gram-Schmidt |
| [`04_minimos_cuadrados.py`](04_minimos_cuadrados.py) | sección 5.4 mínimos cuadrados y ajuste de modelos lineales |
| [`99_verificacion_final.py`](99_verificacion_final.py) | verificación final |

## Fragmentos que imprime el libro

| Archivo | Dónde aparece |
|---|---|
| [`alg-u5-l1-prueba.py`](fragmentos/alg-u5-l1-prueba.py) | Laboratorio del Error 5.1 (prueba con codigo) |
| [`alg-u5-l2-prueba.py`](fragmentos/alg-u5-l2-prueba.py) | Laboratorio del Error 5.2 (prueba con codigo) |
| [`alg-u5-l3-prueba.py`](fragmentos/alg-u5-l3-prueba.py) | Laboratorio del Error 5.3 (prueba con codigo) |
| [`alg-u5-r01-verificacion.py`](fragmentos/alg-u5-r01-verificacion.py) | problema resuelto ALG-U5-01 (verificacion) |
| [`alg-u5-r02-verificacion.py`](fragmentos/alg-u5-r02-verificacion.py) | problema resuelto ALG-U5-02 (verificacion) |
| [`alg-u5-r03-verificacion.py`](fragmentos/alg-u5-r03-verificacion.py) | problema resuelto ALG-U5-03 (verificacion) |
| [`alg-u5-r04-verificacion.py`](fragmentos/alg-u5-r04-verificacion.py) | problema resuelto ALG-U5-04 (verificacion) |
| [`alg-u5-r05-verificacion.py`](fragmentos/alg-u5-r05-verificacion.py) | problema resuelto ALG-U5-05 (verificacion) |
| [`alg-u5-s1-ortogonalidad.py`](fragmentos/alg-u5-s1-ortogonalidad.py) | subtema 5.1: producto interior, angulo y complemento ortogonal en R^n |
| [`alg-u5-s2-proyeccion.py`](fragmentos/alg-u5-s2-proyeccion.py) | subtema 5.2: proyeccion sobre un plano con las ecuaciones normales |
| [`alg-u5-s3-gram-schmidt.py`](fragmentos/alg-u5-s3-gram-schmidt.py) | subtema 5.3: Gram-Schmidt paso a paso y factorizacion QR |
| [`alg-u5-s4-minimos.py`](fragmentos/alg-u5-s4-minimos.py) | subtema 5.4: recta de minimos cuadrados del resorte y sus residuos |

Para ejecutar un fragmento en tu computadora: `python fragmentos/<archivo>.py`. En Colab, copia su contenido en una celda de código.
