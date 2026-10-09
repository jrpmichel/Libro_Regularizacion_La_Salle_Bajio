# ID: INT-U5-R04
# Libro: Int. U5 · problema resuelto INT-U5-04 (verificación)
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-r04-verificacion.py
import sympy as sp

s, t = sp.symbols("s t", nonnegative=True)
tau = 4
integral = sp.integrate(s * sp.exp(s / tau), (s, 0, t))
y = 20 + sp.Rational(3, 2) / tau * sp.exp(-t / tau) * integral
print("y(t) =", sp.expand(sp.simplify(y)))
print("y(20) =", round(float(y.subs(t, 20)), 2))             # 44.04
error = 20 + sp.Rational(3, 2) * t - y
print("error a los 20 s:", round(float(error.subs(t, 20)), 2))   # 5.96
print("error final:", sp.limit(error, t, sp.oo))              # 6
