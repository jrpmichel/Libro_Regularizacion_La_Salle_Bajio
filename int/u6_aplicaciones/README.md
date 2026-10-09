# Int. U6 · Aplicaciones a ingeniería

[![Abrir en Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jrpmichel/Libro_Regularizacion_La_Salle_Bajio/blob/main/int/u6_aplicaciones.ipynb)

Notebook completo: [`u6_aplicaciones.ipynb`](../u6_aplicaciones.ipynb). Cada celda de código del notebook viene de uno de los archivos `NN_*.py` de esta carpeta, en orden. Los archivos de [`fragmentos/`](fragmentos) son los programas que el libro imprime; cada uno corre solo, sin el notebook.

## Archivos del notebook

| Archivo | Contenido |
|---|---|
| [`00_utilidades.py`](00_utilidades.py) | utilidades compartidas |
| [`01_area.py`](01_area.py) | sección 6.1 área entre curvas |
| [`02_volumenes.py`](02_volumenes.py) | sección 6.2 volúmenes de sólidos de revolución |
| [`03_arco.py`](03_arco.py) | sección 6.3 longitud de arco |
| [`04_centroides.py`](04_centroides.py) | sección 6.4 centroides y momentos de inercia de área |
| [`05_masa.py`](05_masa.py) | sección 6.5 momentos de masa |
| [`06_movimiento.py`](06_movimiento.py) | sección 6.6 de la aceleración a la posición |
| [`07_trabajo.py`](07_trabajo.py) | sección 6.7 trabajo, energía y eficiencia |
| [`08_fluidos.py`](08_fluidos.py) | sección 6.8 fuerza hidrostática y caudal |
| [`09_carga.py`](09_carga.py) | sección 6.9 carga, energía y valor eficaz de una señal |
| [`10_costo.py`](10_costo.py) | sección 6.10 costo acumulado |
| [`99_verificacion_final.py`](99_verificacion_final.py) | verificación final |

## Fragmentos que imprime el libro

| Archivo | Dónde aparece |
|---|---|
| [`int-u6-l1-prueba.py`](fragmentos/int-u6-l1-prueba.py) | Laboratorio del Error 6.1 (prueba con código) |
| [`int-u6-l2-prueba.py`](fragmentos/int-u6-l2-prueba.py) | Laboratorio del Error 6.2 (prueba con código) |
| [`int-u6-l3-prueba.py`](fragmentos/int-u6-l3-prueba.py) | Laboratorio del Error 6.3 (prueba con código) |
| [`int-u6-r01-verificacion.py`](fragmentos/int-u6-r01-verificacion.py) | problema resuelto INT-U6-01 (verificación) |
| [`int-u6-r02-verificacion.py`](fragmentos/int-u6-r02-verificacion.py) | problema resuelto INT-U6-02 (verificación) |
| [`int-u6-r03-verificacion.py`](fragmentos/int-u6-r03-verificacion.py) | problema resuelto INT-U6-03 (verificación) |
| [`int-u6-r04-verificacion.py`](fragmentos/int-u6-r04-verificacion.py) | problema resuelto INT-U6-04 (verificación) |
| [`int-u6-r05-verificacion.py`](fragmentos/int-u6-r05-verificacion.py) | problema resuelto INT-U6-05 (verificación) |
| [`int-u6-r06-verificacion.py`](fragmentos/int-u6-r06-verificacion.py) | problema resuelto INT-U6-06 (verificación) |
| [`int-u6-r07-verificacion.py`](fragmentos/int-u6-r07-verificacion.py) | problema resuelto INT-U6-07 (verificación) |
| [`int-u6-r08-verificacion.py`](fragmentos/int-u6-r08-verificacion.py) | problema resuelto INT-U6-08 (verificación) |
| [`int-u6-r09-verificacion.py`](fragmentos/int-u6-r09-verificacion.py) | problema resuelto INT-U6-09 (verificación) |
| [`int-u6-r10-verificacion.py`](fragmentos/int-u6-r10-verificacion.py) | problema resuelto INT-U6-10 (verificación) |
| [`int-u6-s1-area.py`](fragmentos/int-u6-s1-area.py) | subtema 6.1: área entre sen x y cos x en [0, pi] |
| [`int-u6-s10-costo.py`](fragmentos/int-u6-s10-costo.py) | subtema 6.10: costo desde el marginal y curva de aprendizaje |
| [`int-u6-s2-discos.py`](fragmentos/int-u6-s2-discos.py) | subtema 6.2: volumen de un cono sumando discos |
| [`int-u6-s3-poligonal.py`](fragmentos/int-u6-s3-poligonal.py) | subtema 6.3: poligonales sobre y = x^(3/2) |
| [`int-u6-s4-centroide.py`](fragmentos/int-u6-s4-centroide.py) | subtema 6.4: centroide de una región e inercia de un rectángulo |
| [`int-u6-s5-masa.py`](fragmentos/int-u6-s5-masa.py) | subtema 6.5: masa, centro de masa e inercia por trozos |
| [`int-u6-s6-sesgo.py`](fragmentos/int-u6-s6-sesgo.py) | subtema 6.6: doble integración con y sin sesgo |
| [`int-u6-s7-trabajo.py`](fragmentos/int-u6-s7-trabajo.py) | subtema 6.7: trabajo de un resorte y de un bombeo |
| [`int-u6-s8-fluidos.py`](fragmentos/int-u6-s8-fluidos.py) | subtema 6.8: placa sumergida y caudal laminar |
| [`int-u6-s9-carga.py`](fragmentos/int-u6-s9-carga.py) | subtema 6.9: carga total y valor eficaz de un PWM |

Para ejecutar un fragmento en tu computadora: `python fragmentos/<archivo>.py`. En Colab, copia su contenido en una celda de código.
