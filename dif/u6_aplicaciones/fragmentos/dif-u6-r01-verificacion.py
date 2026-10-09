# ID: DIF-U6-R01
# Libro: Dif. U6 · problema resuelto DIF-U6-01 (verificación)
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-r01-verificacion.py
import sympy as sp

t = sp.symbols("t", nonnegative=True)        # min
y1 = 1 - sp.exp(-t)                          # primer orden
y2 = 1 - 2 * sp.exp(-t) + sp.exp(-2 * t)     # segundo orden
print(sp.factor(sp.diff(y2, t)), "|", sp.factor(sp.diff(y2, t, 2)))
ti = sp.solve(sp.diff(y2, t, 2), t)[0]
print("inflexión en t =", ti, "=", round(float(ti), 3), "min")   # ln 2
print("y =", y2.subs(t, ti), "  y' =", sp.diff(y2, t).subs(t, ti))
print("primer orden: y'' =", sp.diff(y1, t, 2), "(siempre negativa)")
