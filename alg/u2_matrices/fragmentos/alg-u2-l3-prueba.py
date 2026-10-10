# ID: ALG-U2-L3
# Libro: Alg. U2 · Laboratorio del Error 2.3 (prueba con codigo)
# Repositorio: alg/u2_matrices/fragmentos/alg-u2-l3-prueba.py
import numpy as np

C = np.array([[2.0, 1.0], [0.0, 1.0]])    # calibracion
R = np.array([[0.0, -1.0], [1.0, 0.0]])   # giro de 90 grados
x = np.array([1.0, 2.0])
y = R @ C @ x                             # primero C, despues R
inv = np.linalg.inv
print("y =", y)                           # [-2.  4.]
print("IA:   R^-1 C^-1 y =", inv(R) @ inv(C) @ y)    # [4. 3.]
print("bien: C^-1 R^-1 y =", inv(C) @ inv(R) @ y)    # [1. 2.]
print("(RC)^-1 == C^-1 R^-1:", np.allclose(inv(R @ C), inv(C) @ inv(R)))
