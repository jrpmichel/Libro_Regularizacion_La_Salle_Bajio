# ID: ALG-U4-L1
# Libro: Alg. U4 · Laboratorio del Error 4.1 (prueba con codigo)
# Repositorio: alg/u4_espacios/fragmentos/alg-u4-l1-prueba.py
import sympy as sp

A = sp.Matrix([[1, 0, 1], [2, 1, 3], [3, 1, 4]])     # capa lineal y = A x
R, piv = A.rref()
print("con la forma reducida:", [list(R[:, j]) for j in piv])   # error
print("con A original:      ", [list(A[:, j]) for j in piv])    # correcto
print("salida con x = (1, 0, 0):", list(A[:, 0]))     # y3 = 3, no 0
print("normal de Col(A):", list(A[:, 0].cross(A[:, 1])))   # y3 = y1 + y2
