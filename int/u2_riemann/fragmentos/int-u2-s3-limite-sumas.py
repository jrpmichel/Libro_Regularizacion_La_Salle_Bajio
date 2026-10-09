# ID: INT-U2-S3
# Libro: Int. U2 · subtema 2.3: la integral de x^2 como límite de sumas derechas
# Repositorio: int/u2_riemann/fragmentos/int-u2-s3-limite-sumas.py
import sympy as sp

k, n = sp.symbols("k n", integer=True, positive=True)
b = sp.Symbol("b", positive=True)
f = lambda x: x**2
dx = b / n
R_n = sp.simplify(sp.summation(f(k * dx) * dx, (k, 1, n)))   # suma derecha
print("R_n =", sp.factor(R_n))
print("límite cuando n -> oo:", sp.limit(R_n, n, sp.oo))     # b**3/3
print("con b = 1 y n = 4:", R_n.subs({b: 1, n: 4}))           # 15/32
