# ID: ALG-U3-S3
# Libro: Alg. U3 · subtema 3.3: Ax como combinacion de columnas y rango
# Repositorio: alg/u3_sistemas/fragmentos/alg-u3-s3-rango.py
import numpy as np
import sympy as sp

A = np.array([[2, 1], [1, 3]])
x = np.array([2, 1])
print("A x =", A @ x)                                 # [5 5]
print("2 a1 + 1 a2 =", 2 * A[:, 0] + 1 * A[:, 1])     # lo mismo

B = sp.Matrix([[1, 2, 3], [2, 4, 7], [3, 6, 10]])
R, pivotes = B.rref()                                 # forma reducida
print("forma reducida:", R.tolist())
print("columnas pivote:", [p + 1 for p in pivotes])   # [1, 3]
print("rango (sympy):", len(pivotes))                 # 2
print("rango (numpy):", np.linalg.matrix_rank(np.array(B, float)))
print("det B =", B.det())                             # 0
