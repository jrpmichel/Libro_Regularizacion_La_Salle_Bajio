# ID: INT-U1-S2
# Libro: Int. U1 · subtema 1.2: sumas inferior y superior de x^2 en [0, 1]
# Repositorio: int/u1_area/fragmentos/int-u1-s2-inferior-superior.py
import numpy as np

f = lambda x: x**2                  # creciente en [0, 1]
a, b = 0, 1
print("    n        s_n        S_n   S_n - s_n")
for n in (4, 10, 100, 1000):
    x = np.linspace(a, b, n + 1)    # n + 1 puntos, n subintervalos
    dx = (b - a) / n
    s = np.sum(f(x[:-1])) * dx      # creciente: mínimo a la izquierda
    S = np.sum(f(x[1:])) * dx       # máximo a la derecha
    print(f"{n:5d}  {s:.6f}  {S:.6f}   {S - s:.6f}")
