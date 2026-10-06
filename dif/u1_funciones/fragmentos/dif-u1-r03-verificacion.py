# ID: DIF-U1-R03
# Libro: Dif. U1 · problema resuelto DIF-U1-03 (verificación)
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-r03-verificacion.py
import sympy as sp
from sympy.calculus.util import continuous_domain

x = sp.symbols("x", real=True)
f = sp.sqrt(x + 3) / (x**2 - 4)
print(continuous_domain(f, x, sp.S.Reals))
print(f.subs(x, -3), f.subs(x, 0))                # 0  -sqrt(3)/4
