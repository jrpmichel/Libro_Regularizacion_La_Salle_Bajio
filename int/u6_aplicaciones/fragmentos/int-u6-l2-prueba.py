# ID: INT-U6-L2
# Libro: Int. U6 · Laboratorio del Error 6.2 (prueba con código)
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-l2-prueba.py
import sympy as sp

x = sp.symbols("x", nonnegative=True)
R, r = sp.sqrt(x), x / 2                     # radios exterior e interior
mal = sp.pi * sp.integrate((R - r)**2, (x, 0, 4))
bien = sp.pi * sp.integrate(R**2 - r**2, (x, 0, 4))
print("(R - r)^2:", mal, "=", round(float(mal), 3))       # 8*pi/15
print("R^2 - r^2:", bien, "=", round(float(bien), 3))     # 8*pi/3
# sólido exterior menos sólido interior
ext = sp.pi * sp.integrate(R**2, (x, 0, 4))
inter = sp.pi * sp.integrate(r**2, (x, 0, 4))
print("exterior - interior:", ext, "-", inter, "=", ext - inter)
