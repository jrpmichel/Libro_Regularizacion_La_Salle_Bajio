# ID: ALG-U3-R01
# Libro: Alg. U3 · problema resuelto ALG-U3-01 (verificacion)
# Repositorio: alg/u3_sistemas/fragmentos/alg-u3-r01-verificacion.py
import sympy as sp

V = [1, 3, 5]                                       # V
h = [sp.Rational(8, 10), sp.Rational(26, 10), sp.Rational(52, 10)]  # m
Ab = sp.Matrix([[1, v, v**2, hv] for v, hv in zip(V, h)])
R, piv = Ab.rref()                                  # Gauss-Jordan exacto
a, b, c = R[:, -1]
print("a, b, c =", a, b, c)                         # 1/5, 1/2, 1/10
for v, hv in zip(V, h):
    print(v, "V ->", a + b * v + c * v**2, "m; medido", hv)
print("h(4) =", float(a + 4 * b + 16 * c), "m")     # 3.8
