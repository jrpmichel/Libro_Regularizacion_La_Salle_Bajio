# ID: DIF-U1-L3
# Libro: Dif. U1 · Laboratorio del Error 1.3 (prueba con código)
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-l3-prueba.py
import sympy as sp
from sympy.calculus.util import continuous_domain

x = sp.symbols("x", real=True)
original = (x**2 - 4) / (x - 2)
print(continuous_domain(original, x, sp.S.Reals))    # sin el 2
print(continuous_domain(x + 2, x, sp.S.Reals))       # todos los reales
print(sp.simplify(original))                         # x + 2
f = lambda v: (v**2 - 4) / (v - 2)
print(round(f(1.999), 3), round(f(2.001), 3))        # 3.999 4.001
try:
    f(2)
except ZeroDivisionError as error:
    print("f(2):", error)                            # la original no existe en 2
