# ID: INT-U4-S2
# Libro: Int. U4 · subtema 4.2: F(b) - F(a) contra el límite de sumas
# Repositorio: int/u4_tfc/fragmentos/int-u4-s2-evaluacion.py
import sympy as sp

x, k, n = sp.symbols("x k n", positive=True)
f = x**2
F = sp.integrate(f, x)                       # x**3/3
print("F(3) - F(0) =", F.subs(x, 3) - F.subs(x, 0))         # 9

dx = 3 / n                                   # sumas derechas de la Unidad 2
Rn = sp.summation(f.subs(x, k * dx) * dx, (k, 1, n))
print("R_n =", sp.factor(Rn))
print("límite de R_n =", sp.limit(Rn, n, sp.oo))            # 9
