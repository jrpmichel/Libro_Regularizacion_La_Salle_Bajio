# ID: ALG-U3-L3
# Libro: Alg. U3 · Laboratorio del Error 3.3 (prueba con codigo)
# Repositorio: alg/u3_sistemas/fragmentos/alg-u3-l3-prueba.py
import numpy as np

A = np.array([[2, 1, 3], [1, 2, 1], [3, 3, 4]], float)
b = np.array([50, 40, 95], float)                   # 95 h de inspeccion
x, *_ = np.linalg.lstsq(A, b, rcond=None)           # no avisa nada
print("lstsq:", x.round(2))                         # [ 7.38 13.19  7.9 ]
print("A x =", (A @ x).round(2), " b =", b)
print("residuo b - A x:", (b - A @ x).round(2))     # no es cero
Ab = np.column_stack([A, b])
print("rangos:", np.linalg.matrix_rank(A), np.linalg.matrix_rank(Ab))
plan = np.array([7, 13, 8])                         # plan redondeado
print("el plan pide:", A @ plan, "h; hay", b)        # 51, 41, 92
