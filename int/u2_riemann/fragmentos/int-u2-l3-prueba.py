# ID: INT-U2-L3
# Libro: Int. U2 · Laboratorio del Error 2.3 (prueba con código)
# Repositorio: int/u2_riemann/fragmentos/int-u2-l3-prueba.py
import sympy as sp

x = sp.symbols("x")
F = 700 - 1000 * x                                 # N
a, b = sp.Rational(4, 10), sp.Rational(1, 10)       # de 0.4 a 0.1 m
print("de 0.1 a 0.4:", sp.integrate(F, (x, b, a)), "J")          # 135
print("de 0.4 a 0.1:", sp.integrate(F, (x, a, b)), "J")          # -135
n = 1000
dx = (0.1 - 0.4) / n                       # negativo: el pistón retrocede
xm = [0.4 + (k + 0.5) * dx for k in range(n)]      # puntos medios
suma = sum((700 - 1000 * xk) * dx for xk in xm)
print("suma en el orden del movimiento:", round(suma, 6), "J")
