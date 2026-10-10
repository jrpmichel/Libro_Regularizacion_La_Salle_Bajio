# ID: ALG-U2-L1
# Libro: Alg. U2 · Laboratorio del Error 2.1 (prueba con codigo)
# Repositorio: alg/u2_matrices/fragmentos/alg-u2-l1-prueba.py
import numpy as np

S = np.array([[2, 0], [0, 1]])      # primero: estira x al doble
R = np.array([[0, -1], [1, 0]])     # despues: gira 90 grados
p = np.array([1, 1])
print("bien, R S p =", R @ S @ p)   # [-1  2]
print("mal,  S R p =", S @ R @ p)   # [-2  1]
print("det(RS) =", round(np.linalg.det(R @ S), 12))   # 2.0
print("det(SR) =", round(np.linalg.det(S @ R), 12))   # 2.0
print("RS == SR:", np.array_equal(R @ S, S @ R))     # False
