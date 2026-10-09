# ID: INT-U4-S5
# Libro: Int. U4 · subtema 4.5: tabla de comprobación contra sumas
# Repositorio: int/u4_tfc/fragmentos/int-u4-s5-comprobacion.py
import numpy as np

f = lambda x: np.cos(2 * x)
a, b = 0.0, 1.0
for nombre, F in (("bien", lambda x: np.sin(2 * x) / 2),
                  ("mal ", lambda x: np.sin(2 * x))):
    valor = F(b) - F(a)
    print(nombre, "F(b) - F(a) =", round(valor, 6))
    for n in (4, 16, 64):
        xs = np.linspace(a, b, n + 1)
        M = f((xs[:-1] + xs[1:]) / 2).sum() * (b - a) / n
        d = M - valor
        print(f"   n = {n:3d}   M_n = {M:.6f}   diferencia = {d:+.2e}")
