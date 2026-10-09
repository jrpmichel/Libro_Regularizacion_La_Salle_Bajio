# ID: INT-U4-R03
# Libro: Int. U4 · problema resuelto INT-U4-03 (verificación)
# Repositorio: int/u4_tfc/fragmentos/int-u4-r03-verificacion.py
import numpy as np
import sympy as sp

t, c = sp.symbols("t c", positive=True)
T = 20 + 180 * (1 - sp.exp(-t / 10))              # grados C, t en min
media = sp.integrate(T, (t, 0, 30)) / 30
print("valor medio:", round(float(media), 2))      # 142.99
L = np.array([float(T.subs(t, v)) for v in (0, 10, 20, 30)])
print("lecturas:", L.round(1))
print("promedio simple:", round(L.mean(), 1))      # 130.1
print("extremos:", round((L[0] + L[-1]) / 2, 1))   # 105.5
print("trapecio:", round((L[0]/2 + L[1] + L[2] + L[3]/2) / 3, 1))  # 138.3
cval = sp.nsolve(T.subs(t, c) - media, c, 10)
print("c =", round(float(cval), 2), "min")         # 11.5
