# ID: DIF-U3-R04
# Libro: Dif. U3 · problema resuelto DIF-U3-04 (verificación)
# Repositorio: dif/u3_transformaciones/fragmentos/dif-u3-r04-verificacion.py
import numpy as np

sig = lambda u: 1 / (1 + np.exp(-u))
neuronas = {"N1": (2, -6), "N2": (0.5, -1.5), "N3": (2, -12)}   # sig(w x + c)
for nombre, (w, c) in neuronas.items():
    h = -c / w                         # umbral: la salida vale 0.5 en x = h
    ancho = 2 * np.log(9) / w          # tramo donde la salida va de 0.1 a 0.9
    print(f"{nombre}: umbral {h}, ancho {ancho:.3f}, salida {sig(w*h + c)}")
print(round(sig(np.log(9)), 3), round(sig(-np.log(9)), 3))      # 0.9 0.1
