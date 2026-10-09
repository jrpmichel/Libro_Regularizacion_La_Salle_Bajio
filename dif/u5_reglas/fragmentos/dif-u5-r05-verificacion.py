# ID: DIF-U5-R05
# Libro: Dif. U5 · problema resuelto DIF-U5-05 (verificación)
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-r05-verificacion.py
import sympy as sp

t = sp.symbols("t", nonnegative=True)       # s; ángulos en radianes
y = 1 - sp.exp(-2 * t) * (sp.cos(4 * t) + sp.sin(4 * t) / 2)
dy = sp.simplify(sp.diff(y, t))
print(dy)                                   # 5*exp(-2*t)*sin(4*t)
print(dy.subs(t, 0), sp.nsolve(dy, t, 0.7)) # 0 y 0.785398 = pi/4
pico = y.subs(t, sp.pi / 4)
print(sp.simplify(pico), round(float(pico), 4))   # 1 + exp(-pi/2), 1.2079
