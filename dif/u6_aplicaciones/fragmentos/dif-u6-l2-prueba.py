# ID: DIF-U6-L2
# Libro: Dif. U6 · Laboratorio del Error 6.2 (prueba con código)
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-l2-prueba.py
import math
import sympy as sp

x = sp.symbols("x")
f, g = 1 - sp.cos(x), x**2 + x
print("en 0:", f.subs(x, 0), "/", g.subs(x, 0))              # 0/0
f1, g1 = sp.diff(f, x), sp.diff(g, x)
print("tras una vez:", f1.subs(x, 0), "/", g1.subs(x, 0))     # 0/1 = 0
print("sympy:", sp.limit(f / g, x, 0))                        # 0
print("numérico:", (1 - math.cos(1e-4)) / (1e-8 + 1e-4))      # casi 0
