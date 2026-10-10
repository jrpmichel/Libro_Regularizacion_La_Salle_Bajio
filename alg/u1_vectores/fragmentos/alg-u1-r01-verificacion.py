# ID: ALG-U1-R01
# Libro: Alg. U1 · problema resuelto ALG-U1-01 (verificacion)
# Repositorio: alg/u1_vectores/fragmentos/alg-u1-r01-verificacion.py
import numpy as np

P = np.array([1, -2, 4])
Q = np.array([3, 4, 1])
PQ, QP = Q - P, P - Q
print("PQ =", PQ, " QP =", QP)          # [2 6 -3]  [-2 -6 3]
print("|PQ| =", np.linalg.norm(PQ))     # 7.0
print("|QP| =", np.linalg.norm(QP))     # 7.0
print("PQ + QP =", PQ + QP)             # vector cero
print("candidatas:", 2 + 6 - 3, np.sqrt(4 + 36 + 9), np.sqrt(4 + 36 - 9))
