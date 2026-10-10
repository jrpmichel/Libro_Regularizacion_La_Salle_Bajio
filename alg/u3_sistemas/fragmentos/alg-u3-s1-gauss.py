# ID: ALG-U3-S1
# Libro: Alg. U3 · subtema 3.1: eliminacion de Gauss-Jordan paso a paso
# Repositorio: alg/u3_sistemas/fragmentos/alg-u3-s1-gauss.py
from fractions import Fraction as F
import numpy as np

M = [[1, 1, 1, 6], [2, 2, 3, 15], [1, 3, 2, 13]]   # matriz aumentada
M = [[F(x) for x in fila] for fila in M]            # fracciones exactas
n = len(M)
for k in range(n):
    if M[k][k] == 0:                                # pivote cero
        j = next(i for i in range(k + 1, n) if M[i][k] != 0)
        M[k], M[j] = M[j], M[k]
        print(f"R{k+1} <-> R{j+1}")
    if M[k][k] != 1:
        print(f"R{k+1} / {M[k][k]}")
        M[k] = [x / M[k][k] for x in M[k]]          # pivote = 1
    for i in range(n):
        if i != k and M[i][k] != 0:                 # ceros arriba y abajo
            m = M[i][k]
            M[i] = [a - m * b for a, b in zip(M[i], M[k])]
            print(f"R{i+1} - ({m}) R{k+1}")
print("solucion:", [str(fila[-1]) for fila in M])  # 1, 2, 3
A = np.array([[1, 1, 1], [2, 2, 3], [1, 3, 2]], float)
print("numpy:   ", np.linalg.solve(A, [6, 15, 13]))
