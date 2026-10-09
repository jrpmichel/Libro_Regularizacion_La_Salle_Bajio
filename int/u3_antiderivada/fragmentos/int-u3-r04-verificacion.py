# ID: INT-U3-R04
# Libro: Int. U3 · problema resuelto INT-U3-04 (verificación)
# Repositorio: int/u3_antiderivada/fragmentos/int-u3-r04-verificacion.py
import sympy as sp

x, C1, C2 = sp.symbols("x C1 C2")
L, P, EI = 3, 4, 2000                     # m, kN, kN m^2
ypp = -P * (L - x) / EI                   # y'' = M/EI
yp = sp.integrate(ypp, x) + C1
yp = yp.subs(C1, sp.solve(yp.subs(x, 0), C1)[0])      # y'(0) = 0
y = sp.integrate(yp, x) + C2
y = y.subs(C2, sp.solve(y.subs(x, 0), C2)[0])         # y(0) = 0
print("y(x) =", sp.factor(y))
print("giro y'(L) =", yp.subs(x, L), "rad")           # -9/1000
print("flecha y(L) =", y.subs(x, L), "m")             # -9/500
print("P L^3/(3 EI) =", sp.Rational(P * L**3, 3 * EI), "m")
