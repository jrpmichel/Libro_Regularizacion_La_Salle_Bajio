# ID: DIF-U6-R07
# Libro: Dif. U6 · problema resuelto DIF-U6-07 (verificación)
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-r07-verificacion.py
import sympy as sp

r = sp.symbols("r", positive=True)               # m
V = sp.Rational(1, 2)                            # m^3
h = V / (sp.pi * r**2)
C = 675 * 2 * sp.pi * r**2 + 450 * 2 * sp.pi * r * h   # pesos
ro = sp.solve(sp.diff(C, r), r)[0]
print("r* =", round(float(ro), 4), "m | h* =", round(float(h.subs(r, ro)), 4),
      "m | h/r =", sp.nsimplify(h.subs(r, ro) / ro))
print("C(r*) =", round(float(C.subs(r, ro)), 2), "pesos")       # 1796.40
print("C'' > 0:", sp.diff(C, r, 2).subs(r, ro) > 0)
print("C(0.35) =", round(float(C.subs(r, 0.35)), 2), "pesos")   # 1805.26
