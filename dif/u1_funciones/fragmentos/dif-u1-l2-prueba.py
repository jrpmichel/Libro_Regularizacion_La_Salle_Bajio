# ID: DIF-U1-L2
# Libro: Dif. U1 · Laboratorio del Error 1.2 (prueba con código)
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-l2-prueba.py
import numpy as np

h = lambda t: 1.5 + 19.6*t - 4.9*t**2                # altura en m, t en s
t_escena = (19.6 + np.sqrt(21.56)) / 9.8
print(round(t_escena, 2), round(h(t_escena), 2))     # 2.47 20.0: la comprobacion pasa
raices = np.sort(np.roots([-4.9, 19.6, 1.5 - 20]))   # h(t) = 20
print(raices.round(2))                               # [1.53 2.47]
print(h(raices).round(2))                            # [20. 20.]
print(raices.sum() / 2)                              # 2.0: eje de simetria de la parabola
