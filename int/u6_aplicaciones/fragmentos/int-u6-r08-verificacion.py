# ID: INT-U6-R08
# Libro: Int. U6 · problema resuelto INT-U6-08 (verificación)
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-r08-verificacion.py
import sympy as sp

h = sp.symbols("h", positive=True)
rho, g, w = 1000, sp.Rational(981, 100), 2       # kg/m^3, m/s^2, m
p = rho * g * h                                   # presión (Pa)
h1, h2 = 1, sp.Rational(5, 2)                     # profundidades (m)
F = sp.integrate(p * w, (h, h1, h2))
hp = sp.integrate(h * p * w, (h, h1, h2)) / F
print("F =", float(F), "N; centro de presión a", round(float(hp), 4), "m")
M = sp.integrate((h - h1) * p * w, (h, h1, h2))   # momento en la bisagra
print("pestillo:", round(float(M / sp.Rational(3, 2))), "N; con el centroide",
      round(float(F * sp.Rational(3, 4) / sp.Rational(3, 2))), "N")
