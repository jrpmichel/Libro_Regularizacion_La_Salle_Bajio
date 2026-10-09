# ID: DIF-U2-L3
# Libro: Dif. U2 · Laboratorio del Error 2.3 (prueba con código)
# Repositorio: dif/u2_limite/fragmentos/dif-u2-l3-prueba.py
import math
import sympy as sp

for x in (10, 1000, 10**5, 10**7):
    print(x, math.sqrt(x**2 + 3*x) - x)            # se acerca a 1.5, no a 0

x = sp.symbols("x", positive=True)
print(sp.limit(sp.sqrt(x**2 + 3*x) - x, x, sp.oo))           # 3/2
racionalizada = 3*x / (sp.sqrt(x**2 + 3*x) + x)
print(sp.limit(racionalizada, x, sp.oo))                      # 3/2
