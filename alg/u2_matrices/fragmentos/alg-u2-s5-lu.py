# ID: ALG-U2-S5
# Libro: Alg. U2 · subtema 2.5: factorizacion LU sin intercambio de renglones
# Repositorio: alg/u2_matrices/fragmentos/alg-u2-s5-lu.py
import numpy as np

def lu(A):
    U = np.array(A, dtype=float)
    n = len(U)
    L = np.eye(n)
    for k in range(n - 1):
        if U[k, k] == 0:
            raise ValueError("pivote cero: hay que intercambiar renglones")
        for i in range(k + 1, n):
            L[i, k] = U[i, k] / U[k, k]          # multiplicador
            U[i] = U[i] - L[i, k] * U[k]         # elimina debajo del pivote
    return L, U

A = [[2, 1, 1], [4, 3, 3], [8, 7, 9]]
L, U = lu(A)
print("L =\n", L, "\nU =\n", U)
print("LU = A:", np.allclose(L @ U, A))
print("det A = producto de la diagonal de U =", np.prod(np.diag(U)))   # 4
