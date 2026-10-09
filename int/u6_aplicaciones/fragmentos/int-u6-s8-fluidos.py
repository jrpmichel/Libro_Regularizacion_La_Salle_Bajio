# ID: INT-U6-S8
# Libro: Int. U6 · subtema 6.8: placa sumergida y caudal laminar
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-s8-fluidos.py
import sympy as sp

h, r, R, vm = sp.symbols("h r R v_max", positive=True)
rho, g = 1000, sp.Rational(981, 100)
# placa cuadrada de 1 m con el borde superior en la superficie
F = sp.integrate(rho * g * h * 1, (h, 0, 1))
hp = sp.integrate(h * rho * g * h, (h, 0, 1)) / F
print("F =", F, "N; centro de presión a", hp, "m")       # 4905; 2/3
# tubo de radio R con perfil laminar
Q = sp.integrate(vm * (1 - r**2 / R**2) * 2 * sp.pi * r, (r, 0, R))
print("Q =", sp.simplify(Q), "; media:", sp.simplify(Q / (sp.pi * R**2)))
