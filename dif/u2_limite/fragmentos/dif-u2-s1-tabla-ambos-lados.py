# ID: DIF-U2-S1
# Libro: Dif. U2 · subtema 2.1: tabla desde ambos lados
# Repositorio: dif/u2_limite/fragmentos/dif-u2-s1-tabla-ambos-lados.py
import numpy as np
import matplotlib.pyplot as plt

f = lambda x: (x**3 - 1) / (x - 1)          # no está definida en x = 1
for h in (0.1, 0.01, 0.001):
    print(f"x = {1-h:<6} f = {f(1-h):.6f}  |  x = {1+h:<6} f = {f(1+h):.6f}")

xx = np.linspace(0, 2, 400)
xx = xx[np.abs(xx - 1) > 1e-9]                # se excluye x = 1
xl = np.array([0.5, 0.75, 0.9]); xr = np.array([1.5, 1.25, 1.1])
plt.plot(xx, f(xx), color="gray")
plt.plot(xl, f(xl), "^", label="por la izquierda")
plt.plot(xr, f(xr), "s", mfc="white", label="por la derecha")
plt.plot(1, 3, "o", mfc="white", mec="red")  # hueco en (1, 3)
plt.xlabel("x"); plt.ylabel("y"); plt.legend(); plt.grid(True); plt.show()
