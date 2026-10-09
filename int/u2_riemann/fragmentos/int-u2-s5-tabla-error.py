# ID: INT-U2-S5
# Libro: Int. U2 · subtema 2.5: tabla de error contra n de cuatro reglas
# Repositorio: int/u2_riemann/fragmentos/int-u2-s5-tabla-error.py
import numpy as np
from math import erf, sqrt, pi

phi = lambda x: np.exp(-x**2 / 2) / np.sqrt(2 * pi)
exacto = 0.5 * erf(1 / sqrt(2))                # P(0 <= Z <= 1) = 0.3413447...

def reglas(f, a, b, n):
    x = np.linspace(a, b, n + 1); y = f(x); h = (b - a) / n
    L = h * y[:-1].sum()
    M = h * f((x[:-1] + x[1:]) / 2).sum()
    T = h * (y.sum() - (y[0] + y[-1]) / 2)
    S = h / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum())
    return L, M, T, S

print("  n   error L   error M   error T   error S")
for n in (2, 4, 8, 16, 32):
    errores = [abs(v - exacto) for v in reglas(phi, 0, 1, n)]
    print(f"{n:3d}  " + "  ".join(f"{e:.2e}" for e in errores))
