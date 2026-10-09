# ID: DIF-U3-S1
# Libro: Dif. U3 · subtema 3.1: correspondencia de puntos en y = a f(b(x - h)) + k
# Repositorio: dif/u3_transformaciones/fragmentos/dif-u3-s1-correspondencia.py
import numpy as np
import matplotlib.pyplot as plt

f = np.sqrt                                   # función base
a, b, h, k = -2, 2, 1, 3                      # g(x) = a f(b(x - h)) + k
g = lambda x: a * f(b * (x - h)) + k

for x0, y0 in [(0, 0), (1, 1), (4, 2)]:       # puntos de la gráfica de f
    x1, y1 = h + x0 / b, k + a * y0           # su imagen en la gráfica de g
    print(f"({x0}, {y0}) -> ({x1}, {y1})   comprobación: g({x1}) = {g(x1)}")

x = np.linspace(0, 5, 300); xg = np.linspace(1, 5, 300)
plt.plot(x, f(x), label=r"$f(x)=\sqrt{x}$")
plt.plot(xg, g(xg), "--", label=r"$g(x)=-2\sqrt{2(x-1)}+3$")
plt.grid(True); plt.legend(); plt.xlabel("x"); plt.show()
