# ID: ALG-U4-S1
# Libro: Alg. U4 · subtema 4.1: ¿esta b en el generado de dos vectores?
# Repositorio: alg/u4_espacios/fragmentos/alg-u4-s1-generado.py
import numpy as np
import sympy as sp

u, v = sp.Matrix([1, 0, 1]), sp.Matrix([2, 1, 1])
A = u.row_join(v)                                # los vectores como columnas
for b in ([5, 4, 1], [5, 4, 3]):
    b = sp.Matrix(b)
    r, r_aum = A.rank(), A.row_join(b).rank()
    print("b =", list(b), " rangos:", r, r_aum,
          "-> en gen{u, v}" if r == r_aum else "-> fuera")
n = u.cross(v)                                   # normal del plano gen{u, v}
print("normal:", list(n))                        # [-1, 1, 1]
c = np.random.default_rng(1).normal(size=(5, 2))  # 5 combinaciones al azar
puntos = c @ np.array([[1, 0, 1], [2, 1, 1]])
print("n . (c1 u + c2 v):", np.round(puntos @ np.array([-1, 1, 1]), 12))
