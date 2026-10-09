# ID: DIF-U4-S6
# Libro: Dif. U4 · subtema 4.6: cocientes por la izquierda y la derecha
# Repositorio: dif/u4_derivada/fragmentos/dif-u4-s6-laterales.py
import numpy as np

casos = {"|x|": np.abs, "x^(2/3)": lambda x: np.cbrt(x**2),
         "x^(1/3)": np.cbrt}
for nombre, f in casos.items():
    for h in [0.1, 0.001]:
        izq = (f(-h) - f(0)) / (-h)      # incremento -h < 0: izquierda
        der = (f(h) - f(0)) / h          # incremento h > 0: derecha
        print(f"{nombre:8} h={h:<6} izq={izq:9.3f}  der={der:9.3f}")
