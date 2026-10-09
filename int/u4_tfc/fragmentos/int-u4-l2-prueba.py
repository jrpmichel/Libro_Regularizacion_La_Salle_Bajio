# ID: INT-U4-L2
# Libro: Int. U4 · Laboratorio del Error 4.2 (prueba con código)
# Repositorio: int/u4_tfc/fragmentos/int-u4-l2-prueba.py
import sympy as sp

x = sp.symbols("x", real=True)
h = sp.Rational(1, 100) * x**2 - sp.Rational(1, 5) * x   # m
b = 10                                                   # ancho, m
neto = b * sp.integrate(h, (x, 0, 40))
relleno = b * abs(sp.integrate(h, (x, 0, 20)))
corte = b * sp.integrate(h, (x, 20, 40))
print("neto:", float(neto), "m^3")                     # 533.3
print("relleno:", float(relleno), " corte:", float(corte))
print("tierra movida:", float(relleno + corte), "m^3")   # 800
