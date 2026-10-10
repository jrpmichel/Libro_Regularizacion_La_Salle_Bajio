# ID: ALG-U2-L2
# Libro: Alg. U2 · Laboratorio del Error 2.2 (prueba con codigo)
# Repositorio: alg/u2_matrices/fragmentos/alg-u2-l2-prueba.py
import numpy as np

T = 2 * np.eye(3)                   # escala al doble en x, y y z
V0 = 150.0                          # cm3
print("det(2I) =", round(np.linalg.det(T), 10))      # 8, no 2
print("volumen nuevo =", round(np.linalg.det(T) * V0, 6), "cm3")   # 1200
A = np.array([[1.0, 2, 0], [0, 1, 1], [2, 0, 1]])
print("det A =", round(np.linalg.det(A), 10))          # 5.0
print("det 2A =", round(np.linalg.det(2 * A), 10))     # 40.0 = 2**3 * 5
