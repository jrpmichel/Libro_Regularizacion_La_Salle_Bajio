# ID: DIF-U5-R03
# Libro: Dif. U5 · problema resuelto DIF-U5-03 (verificación)
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-r03-verificacion.py
import sympy as sp

t = sp.symbols("t", nonnegative=True)
d = sp.sqrt((12 * t)**2 + 50**2)                 # m
dp = sp.simplify(sp.diff(d, t))
print(dp)                       # 72*t/sqrt(36*t**2 + 625), la misma fracción
print(d.subs(t, 10), dp.subs(t, 10))             # 130 144/13
print(round(float(dp.subs(t, 10)), 3))           # 11.077 m/s
print(dp.subs(t, 0), sp.limit(dp, t, sp.oo))     # 0 y 12
