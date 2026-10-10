# ID: ALG-U3-S2
# Libro: Alg. U3 · subtema 3.2: regla de Cramer con determinantes
# Repositorio: alg/u3_sistemas/fragmentos/alg-u3-s2-cramer.py
import numpy as np

def cramer(A, b):
    A, b = np.array(A, float), np.array(b, float)
    if np.linalg.matrix_rank(A) < len(A):           # det A = 0
        raise ValueError("det A = 0: Cramer no se aplica; usa el rango")
    d = np.linalg.det(A)
    x = []
    for i in range(len(A)):
        Ai = A.copy()
        Ai[:, i] = b                                 # columna i por b
        x.append(np.linalg.det(Ai) / d)
    return np.array(x)

A, b = [[2, 1], [1, 3]], [5, 5]                     # sistema de la figura
print("Cramer:", cramer(A, b).round(12))            # [2. 1.]
print("solve: ", np.linalg.solve(A, b))
try:
    cramer([[1, 2], [2, 4]], [3, 7])                # rectas paralelas
except ValueError as e:
    print(e)
