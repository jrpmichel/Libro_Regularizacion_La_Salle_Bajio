# ID: INT-U2-R04
# Libro: Int. U2 · problema resuelto INT-U2-04 (verificación)
# Repositorio: int/u2_riemann/fragmentos/int-u2-r04-verificacion.py
import numpy as np
from math import erf, sqrt, pi

phi = lambda x: np.exp(-x**2 / 2) / np.sqrt(2 * pi)
P = 0.5 * erf(1 / sqrt(2))                     # 0.3413447...

def T(n):
    y = phi(np.linspace(0, 1, n + 1))
    return (y.sum() - (y[0] + y[-1]) / 2) / n

def S(n):
    y = phi(np.linspace(0, 1, n + 1))
    impares, pares = y[1:-1:2].sum(), y[2:-1:2].sum()
    return (y[0] + y[-1] + 4 * impares + 2 * pares) / (3 * n)

for n in (2, 4, 8, 16):
    print(n, f"{T(n) - P:.3e}", f"{S(n) - P:.3e}")
nT = next(n for n in range(2, 1000) if abs(T(n) - P) < 5e-7)
nS = next(n for n in range(2, 100, 2) if abs(S(n) - P) < 5e-7)
print("trapecio: n =", nT, "   Simpson: n =", nS)       # 201 y 10
print("error de S(10):", f"{S(10) - P:.2e}")
