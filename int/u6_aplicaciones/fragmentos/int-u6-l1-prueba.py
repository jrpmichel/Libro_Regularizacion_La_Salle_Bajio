# ID: INT-U6-L1
# Libro: Int. U6 · Laboratorio del Error 6.1 (prueba con código)
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-l1-prueba.py
import sympy as sp

x = sp.symbols("x", real=True)
f, g = 4 * x, x**3
print("con signo:", sp.integrate(f - g, (x, 0, 3)))          # -9/4
print("cruces:", sp.solve(sp.Eq(f, g), x))                     # -2, 0, 2
area = sp.integrate(f - g, (x, 0, 2)) + sp.integrate(g - f, (x, 2, 3))
print("área:", area, "=", float(area))                         # 41/4
