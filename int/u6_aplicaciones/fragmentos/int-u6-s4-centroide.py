# ID: INT-U6-S4
# Libro: Int. U6 · subtema 6.4: centroide de una región e inercia de un rectángulo
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-s4-centroide.py
import sympy as sp

x, y, b, h = sp.symbols("x y b h", positive=True)
f = 4 - x**2                               # región bajo f en [0, 2]
A = sp.integrate(f, (x, 0, 2))
xc = sp.integrate(x * f, (x, 0, 2)) / A
yc = sp.integrate(f**2 / 2, (x, 0, 2)) / A
print("A =", A, "; centroide:", (xc, yc))  # 16/3; (3/4, 8/5)
# rectángulo b x h con franjas horizontales de ancho b
Ic = sp.integrate(y**2 * b, (y, -h / 2, h / 2))
Ib = sp.integrate(y**2 * b, (y, 0, h))
print("I centroidal:", Ic, "; I en la base:", Ib)  # b h^3/12, b h^3/3
