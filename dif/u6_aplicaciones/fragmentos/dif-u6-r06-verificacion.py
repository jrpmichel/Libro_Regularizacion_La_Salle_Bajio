# ID: DIF-U6-R06
# Libro: Dif. U6 · problema resuelto DIF-U6-06 (verificación)
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-r06-verificacion.py
import sympy as sp

t = sp.symbols("t", nonnegative=True)
D, T = sp.Rational(2, 5), 2                       # m, s
tau = t / T
s = D * (10 * tau**3 - 15 * tau**4 + 6 * tau**5)
v, a = sp.diff(s, t), sp.diff(s, t, 2)
print(sp.factor(v), "|", sp.factor(a))
print("extremos:", s.subs(t, 0), s.subs(t, T), v.subs(t, 0),
      v.subs(t, T), a.subs(t, 0), a.subs(t, T))
print("v máx:", v.subs(t, 1), "m/s")              # 3/8 en t = 1 s
t1 = min(sp.solve(sp.diff(a, t), t))
print("a máx en t =", round(float(t1), 4), ":",
      round(float(a.subs(t, t1)), 4), "m/s^2")    # 0.4226 s, 0.5774
