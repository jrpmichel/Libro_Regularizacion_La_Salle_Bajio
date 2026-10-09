# ID: DIF-U3-L2
# Libro: Dif. U3 · Laboratorio del Error 3.2 (prueba con código)
# Repositorio: dif/u3_transformaciones/fragmentos/dif-u3-l2-prueba.py
import math
import sympy as sp
from sympy.calculus.util import continuous_domain

x = sp.symbols("x", real=True)
print(continuous_domain(sp.sqrt(2*x - 6), x, sp.S.Reals))   # Interval(3, oo)
print(math.sqrt(2*4 - 6))       # x = 4, que es menor que 6, sí está: sqrt(2)
print(sp.factor(2*x - 6))       # 2*(x - 3): el desplazamiento es 3
