# ID: ALG-U5-S1
# Libro: Alg. U5 · subtema 5.1: producto interior, angulo y complemento ortogonal en R^n
# Repositorio: alg/u5_ortogonalidad/fragmentos/alg-u5-s1-ortogonalidad.py
import numpy as np
import sympy as sp

u, v = np.array([1, 2, 2]), np.array([2, 1, -2])
print("u . v =", u @ v)                                   # 0: ortogonales
print("|u|^2 + |v|^2 =", u @ u + v @ v, "  |u + v|^2 =", (u + v) @ (u + v))
a, b = np.array([1, 1, 1, 1]), np.array([1, 0, 0, 0])
cos = a @ b / (np.linalg.norm(a) * np.linalg.norm(b))
print("angulo en R^4:", round(np.degrees(np.arccos(cos)), 6), "grados")   # 60
W = sp.Matrix([[1, 1, 1]])                 # renglon: el vector que genera W
Wperp = W.nullspace()                      # vectores perpendiculares a (1,1,1)
print("base de W^perp:", [list(w) for w in Wperp])
print("dim W + dim W^perp =", W.rank(), "+", len(Wperp), "=", W.cols)
