# ID: INT-U6-S7
# Libro: Int. U6 · subtema 6.7: trabajo de un resorte y de un bombeo
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-s7-trabajo.py
import sympy as sp

x, y = sp.symbols("x y", real=True)
W_res = sp.integrate(800 * x, (x, sp.Rational(5, 100), sp.Rational(15, 100)))
print("resorte de 800 N/m, de 5 a 15 cm:", W_res, "J")       # 8
rho, g, r, H = 1000, 9.81, 1, 3
# rebanada a la altura y: pesa rho g pi r^2 dy y sube H - y
W_bomba = sp.integrate(rho * g * sp.pi * r**2 * (H - y), (y, 0, H))
print("vaciar el tanque:", round(float(W_bomba), 1), "J")     # 138685.6
