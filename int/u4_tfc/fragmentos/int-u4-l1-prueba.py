# ID: INT-U4-L1
# Libro: Int. U4 · Laboratorio del Error 4.1 (prueba con código)
# Repositorio: int/u4_tfc/fragmentos/int-u4-l1-prueba.py
import sympy as sp

x, t = sp.symbols("x t", real=True)
G = sp.integrate(sp.cos(t), (t, 0, x**2))      # sen(x**2)
bien = sp.diff(G, x)
mal = sp.cos(x**2)
print("G(x) =", G, "   G'(x) =", bien)
print("en x = 1: mal", round(float(mal.subs(x, 1)), 3),
      " bien", round(float(bien.subs(x, 1)), 3))   # 0.540 y 1.081
