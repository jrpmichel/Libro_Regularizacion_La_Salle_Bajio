# ID: DIF-U5-R04
# Libro: Dif. U5 · problema resuelto DIF-U5-04 (verificación)
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-r04-verificacion.py
import sympy as sp

x = sp.symbols("x")
w, L, EI = 2, 3, 2000                       # kN/m, m, kN*m^2
y = -sp.Rational(w, 24 * EI) * (x**4 - 4 * L * x**3 + 6 * L**2 * x**2)
d1, d2, d3, d4 = [sp.diff(y, x, n) for n in (1, 2, 3, 4)]
print(float(y.subs(x, L)), float(d1.subs(x, L)))  # -0.010125 m, -0.0045
print(EI * d2.subs(x, 0), EI * d3.subs(x, 0))     # M(0) = -9, V(0) = 6
print(EI * d4)                                    # -2 = -w: la carga
print(EI * d2.subs(x, L), EI * d3.subs(x, L))     # 0 y 0 en x = L
