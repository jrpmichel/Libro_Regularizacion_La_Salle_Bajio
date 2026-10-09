# ID: INT-U2-L1
# Libro: Int. U2 · Laboratorio del Error 2.1 (prueba con código)
# Repositorio: int/u2_riemann/fragmentos/int-u2-l1-prueba.py
import sympy as sp

k, n = sp.symbols("k n", integer=True, positive=True)
print("suma de 1, de k = 1 a n:", sp.summation(1, (k, 1, n)))   # n, no 1
R_n = sp.summation(((k / n)**2 + 1) / n, (k, 1, n))
print("límite correcto:", sp.limit(R_n, n, sp.oo))              # 4/3
print("R_1000 =", float(R_n.subs(n, 1000)))                      # 1.3338...
