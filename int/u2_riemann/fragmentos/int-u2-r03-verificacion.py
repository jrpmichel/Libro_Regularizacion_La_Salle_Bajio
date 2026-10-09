# ID: INT-U2-R03
# Libro: Int. U2 · problema resuelto INT-U2-03 (verificación)
# Repositorio: int/u2_riemann/fragmentos/int-u2-r03-verificacion.py
import sympy as sp

x = sp.symbols("x", real=True)
fondo = x**2 / 2                                    # m
bordes = sp.solve(sp.Eq(fondo, 2), x)               # [-2, 2]
A = sp.integrate(2 - fondo, (x, bordes[0], bordes[1]))
print("bordes:", bordes, "  A =", A, "=", round(float(A), 4), "m^2")
print("Q =", float(sp.Rational(6, 10) * A), "m^3/s")          # 3.2
triangulo = sp.Rational(1, 2) * 4 * 2                         # 4 m^2
print("Arquímedes, 4/3 del triángulo:", sp.Rational(4, 3) * triangulo)
