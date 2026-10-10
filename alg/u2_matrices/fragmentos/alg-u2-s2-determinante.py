# ID: ALG-U2-S2
# Libro: Alg. U2 · subtema 2.2: determinante por Sarrus y por cofactores
# Repositorio: alg/u2_matrices/fragmentos/alg-u2-s2-determinante.py
import numpy as np

def sarrus(A):                      # solo 3x3
    (a, b, c), (d, e, f), (g, h, i) = A
    return a*e*i + b*f*g + c*d*h - c*e*g - a*f*h - b*d*i

def cofactores(A):                  # cualquier n x n, por el renglon 1
    A = np.asarray(A)
    if len(A) == 1:
        return A[0, 0]
    total = 0
    for j in range(len(A)):
        menor = np.delete(np.delete(A, 0, axis=0), j, axis=1)
        total += (-1) ** j * A[0, j] * cofactores(menor)
    return total

A = [[2, -1, 3], [0, 4, 1], [5, 2, -2]]
print("Sarrus:", sarrus(A), " cofactores:", cofactores(A))       # -85 y -85
print("numpy:", round(np.linalg.det(A), 10))                      # -85.0
