# ID: INT-U2-R01
# Libro: Int. U2 · problema resuelto INT-U2-01 (verificación)
# Repositorio: int/u2_riemann/fragmentos/int-u2-r01-verificacion.py
import numpy as np
import sympy as sp

Q = lambda t: 30 / (1 + t)                     # m^3/s, t en horas
t = np.arange(5.0)                             # 0, 1, 2, 3, 4 h; dt = 1 h
L, R, M = Q(t[:-1]).sum(), Q(t[1:]).sum(), Q(t[:-1] + 0.5).sum()
for nombre, v in (("L4", L), ("R4", R), ("M4", M)):
    print(nombre, round(v, 4), "->", round(v * 3600), "m^3")
x = sp.symbols("x")
q = 30 / (1 + x)
print("Q' =", sp.diff(q, x), "   Q'' =", sp.diff(q, x, 2))
exacto = sp.integrate(q, (x, 0, 4))            # 30 log(5)
print("exacto:", exacto, "->", round(float(exacto) * 3600), "m^3")
