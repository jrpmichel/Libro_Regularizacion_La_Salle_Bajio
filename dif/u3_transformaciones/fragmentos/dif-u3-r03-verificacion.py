# ID: DIF-U3-R03
# Libro: Dif. U3 · problema resuelto DIF-U3-03 (verificación)
# Repositorio: dif/u3_transformaciones/fragmentos/dif-u3-r03-verificacion.py
import numpy as np
import sympy as sp

t = sp.symbols("t")
v = 10 * sp.sin(120 * sp.pi * t - sp.pi / 3)
h = sp.Rational(1, 360)                                  # (pi/3)/(120 pi)
print(sp.simplify(v - 10 * sp.sin(120 * sp.pi * (t - h))))   # 0
print("h =", h, "s =", round(1000 * float(h), 3), "ms")      # 2.778
print("primer máximo en t =", h + sp.Rational(1, 240), "s")  # 1/144
tt = np.linspace(0, 1 / 60, 100001)
vv = 10 * np.sin(120 * np.pi * tt - np.pi / 3)
print("máximo numérico:", round(1000 * tt[np.argmax(vv)], 2), "ms")  # 6.94
