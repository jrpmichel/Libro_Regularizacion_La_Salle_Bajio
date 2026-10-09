# ID: DIF-U6-S6
# Libro: Dif. U6 · subtema 6.6: velocidad y aceleración de un auto que frena
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-s6-frenado.py
import sympy as sp

t = sp.symbols("t", nonnegative=True)
s = 20 * t - sp.Rational(5, 2) * t**2       # m, con t en s
v, a = sp.diff(s, t), sp.diff(s, t, 2)
print("v(t) =", v, "m/s   a(t) =", a, "m/s^2")
alto = sp.solve(v, t)[0]                     # se detiene cuando v = 0
print("se detiene en t =", alto, "s, tras", s.subs(t, alto), "m")
