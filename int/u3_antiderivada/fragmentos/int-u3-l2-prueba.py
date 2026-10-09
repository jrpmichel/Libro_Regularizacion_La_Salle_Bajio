# ID: INT-U3-L2
# Libro: Int. U3 · Laboratorio del Error 3.2 (prueba con código)
# Repositorio: int/u3_antiderivada/fragmentos/int-u3-l2-prueba.py
import math
import sympy as sp

t = sp.symbols("t")
L = sp.Rational(1, 2)                         # H (supuesto)
v = 170 * sp.sin(120 * sp.pi * t)             # V
mal = -(170 / L) * sp.cos(120 * sp.pi * t)
bien = sp.integrate(v, t) / L                 # sympy divide entre 120 pi
print("corriente bien:", bien)
print("pico mal:", float(170 / L), "A")
print("pico bien:", round(float(170 / L) / (120 * math.pi), 3), "A")
razon = sp.simplify(sp.diff(mal, t) / (v / L))
print("derivada de la mala entre v/L:", razon)        # 120 pi
