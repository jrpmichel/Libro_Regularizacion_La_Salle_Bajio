# Dif. U4 · La derivada

[![Abrir en Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jrpmichel/Libro_Regularizacion_La_Salle_Bajio/blob/main/dif/u4_derivada.ipynb)

Notebook completo: [`u4_derivada.ipynb`](../u4_derivada.ipynb). Cada celda de código del notebook viene de uno de los archivos `NN_*.py` de esta carpeta, en orden. Los archivos de [`fragmentos/`](fragmentos) son los programas que el libro imprime; cada uno corre solo, sin el notebook.

## Archivos del notebook

| Archivo | Contenido |
|---|---|
| [`00_utilidades.py`](00_utilidades.py) | utilidades compartidas |
| [`01_pendiente_secante.py`](01_pendiente_secante.py) | sección 4.1 pendiente de una recta y de la secante |
| [`02_paso_al_limite.py`](02_paso_al_limite.py) | sección 4.2 paso al límite |
| [`03_definicion.py`](03_definicion.py) | sección 4.3 primeros ejemplos |
| [`04_recta_tangente.py`](04_recta_tangente.py) | sección 4.4 notaciones y recta tangente |
| [`05_razon_cambio.py`](05_razon_cambio.py) | sección 4.5 la derivada como razón de cambio |
| [`06_no_derivables.py`](06_no_derivables.py) | sección 4.6 derivable implica continua; puntos no derivables |
| [`99_verificacion_final.py`](99_verificacion_final.py) | verificación final |

## Fragmentos que imprime el libro

| Archivo | Dónde aparece |
|---|---|
| [`dif-u4-l1-prueba.py`](fragmentos/dif-u4-l1-prueba.py) | Laboratorio del Error 4.1 (prueba con código) |
| [`dif-u4-l2-prueba.py`](fragmentos/dif-u4-l2-prueba.py) | Laboratorio del Error 4.2 (prueba con código) |
| [`dif-u4-l3-prueba.py`](fragmentos/dif-u4-l3-prueba.py) | Laboratorio del Error 4.3 (prueba con código) |
| [`dif-u4-r01-verificacion.py`](fragmentos/dif-u4-r01-verificacion.py) | problema resuelto DIF-U4-01 (verificación) |
| [`dif-u4-r02-verificacion.py`](fragmentos/dif-u4-r02-verificacion.py) | problema resuelto DIF-U4-02 (verificación) |
| [`dif-u4-r03-verificacion.py`](fragmentos/dif-u4-r03-verificacion.py) | problema resuelto DIF-U4-03 (verificación) |
| [`dif-u4-r04-verificacion.py`](fragmentos/dif-u4-r04-verificacion.py) | problema resuelto DIF-U4-04 (verificación) |
| [`dif-u4-s1-pendientes.py`](fragmentos/dif-u4-s1-pendientes.py) | subtema 4.1: pendiente de una recta y de una secante |
| [`dif-u4-s2-paso-al-limite.py`](fragmentos/dif-u4-s2-paso-al-limite.py) | subtema 4.2: secantes del café con h cada vez menor |
| [`dif-u4-s3-definicion.py`](fragmentos/dif-u4-s3-definicion.py) | subtema 4.3: derivada por definición con sympy |
| [`dif-u4-s4-recta-tangente.py`](fragmentos/dif-u4-s4-recta-tangente.py) | subtema 4.4: recta tangente a 1/x en x = 2 |
| [`dif-u4-s5-razon-cambio.py`](fragmentos/dif-u4-s5-razon-cambio.py) | subtema 4.5: velocidad media e instantánea |
| [`dif-u4-s6-laterales.py`](fragmentos/dif-u4-s6-laterales.py) | subtema 4.6: cocientes por la izquierda y la derecha |

Para ejecutar un fragmento en tu computadora: `python fragmentos/<archivo>.py`. En Colab, copia su contenido en una celda de código.
