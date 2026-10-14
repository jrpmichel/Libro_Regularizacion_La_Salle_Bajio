# Álg. U6 · Aplicaciones a ingeniería

[![Abrir en Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jrpmichel/Libro_Regularizacion_La_Salle_Bajio/blob/main/alg/u6_aplicaciones.ipynb)

Notebook completo: [`u6_aplicaciones.ipynb`](../u6_aplicaciones.ipynb). Cada celda de código del notebook viene de uno de los archivos `NN_*.py` de esta carpeta, en orden. Los archivos de [`fragmentos/`](fragmentos) son los programas que el libro imprime; cada uno corre solo, sin el notebook.

## Archivos del notebook

| Archivo | Contenido |
|---|---|
| [`00_utilidades.py`](00_utilidades.py) | utilidades compartidas |
| [`01_equilibrio.py`](01_equilibrio.py) | sección 6.1 equilibrio de fuerzas en 3D, momento y trabajo |
| [`02_mallas.py`](02_mallas.py) | sección 6.2 análisis de mallas con leyes de Kirchhoff |
| [`03_mezclas.py`](03_mezclas.py) | sección 6.3 mezclas y asignación de recursos |
| [`04_armaduras.py`](04_armaduras.py) | sección 6.4 armadura por el método de nodos |
| [`05_ajuste.py`](05_ajuste.py) | sección 6.5 ajuste de datos por mínimos cuadrados |
| [`06_transformaciones.py`](06_transformaciones.py) | sección 6.6 transformaciones de coordenadas |
| [`99_verificacion_final.py`](99_verificacion_final.py) | verificación final |

## Fragmentos que imprime el libro

| Archivo | Dónde aparece |
|---|---|
| [`alg-u6-l1-prueba.py`](fragmentos/alg-u6-l1-prueba.py) | Laboratorio del Error: el momento al reves |
| [`alg-u6-l2-prueba.py`](fragmentos/alg-u6-l2-prueba.py) | Laboratorio del Error: el balance que siempre cierra |
| [`alg-u6-l3-prueba.py`](fragmentos/alg-u6-l3-prueba.py) | Laboratorio del Error: la inversa que solo probo traslaciones |
| [`alg-u6-r01-verificacion.py`](fragmentos/alg-u6-r01-verificacion.py) | problema resuelto ALG-U6-01 (verificacion) |
| [`alg-u6-r02-verificacion.py`](fragmentos/alg-u6-r02-verificacion.py) | problema resuelto ALG-U6-02 (verificacion) |
| [`alg-u6-r03-verificacion.py`](fragmentos/alg-u6-r03-verificacion.py) | problema resuelto ALG-U6-03 (verificacion) |
| [`alg-u6-r04-verificacion.py`](fragmentos/alg-u6-r04-verificacion.py) | problema resuelto ALG-U6-04 (verificacion) |
| [`alg-u6-r05-verificacion.py`](fragmentos/alg-u6-r05-verificacion.py) | problema resuelto ALG-U6-05 (verificacion) |
| [`alg-u6-r06-verificacion.py`](fragmentos/alg-u6-r06-verificacion.py) | problema resuelto ALG-U6-06 (verificacion) |
| [`alg-u6-r07-verificacion.py`](fragmentos/alg-u6-r07-verificacion.py) | problema resuelto ALG-U6-07 (verificacion) |
| [`alg-u6-s1-equilibrio.py`](fragmentos/alg-u6-s1-equilibrio.py) | subtema 6.1: tensiones de tres cables, momento y trabajo |
| [`alg-u6-s2-mallas.py`](fragmentos/alg-u6-s2-mallas.py) | subtema 6.2: dos mallas, signos de las fuentes y comprobacion |
| [`alg-u6-s3-mezclas.py`](fragmentos/alg-u6-s3-mezclas.py) | subtema 6.3: mezclas posibles con tres chatarras |
| [`alg-u6-s4-armadura.py`](fragmentos/alg-u6-s4-armadura.py) | subtema 6.4: metodo de nodos para la armadura de la semilla |
| [`alg-u6-s5-ajuste.py`](fragmentos/alg-u6-s5-ajuste.py) | subtema 6.5: recta por el origen y modelo de primer orden |
| [`alg-u6-s6-transformaciones.py`](fragmentos/alg-u6-s6-transformaciones.py) | subtema 6.6: matrices homogeneas, inversa y orden |

Para ejecutar un fragmento en tu computadora: `python fragmentos/<archivo>.py`. En Colab, copia su contenido en una celda de código.
