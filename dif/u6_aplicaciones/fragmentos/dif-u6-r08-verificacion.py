# ID: DIF-U6-R08
# Libro: Dif. U6 · problema resuelto DIF-U6-08 (verificación)
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-r08-verificacion.py
import sympy as sp

t = sp.symbols("t", nonnegative=True)                 # s
v = 50 * (sp.exp(-500 * t) - sp.exp(-2000 * t))       # V
i = 2e-6 * sp.diff(v, t)                              # A
tp = sp.log(4) / 1500                                 # donde i = 0
print("v'(tp) =", sp.simplify(sp.diff(v, t).subs(t, tp)))
print("tp =", round(float(tp) * 1e3, 3), "ms | v máx =",
      round(float(v.subs(t, tp)), 2), "V")            # 0.924 ms, 23.62 V
print("i(0) =", float(i.subs(t, 0)), "A")             # 0.15 A
print("i mín en", round(2 * float(tp) * 1e3, 3), "ms:",
      round(float(i.subs(t, 2 * tp)), 4), "A")         # -0.0149 A
