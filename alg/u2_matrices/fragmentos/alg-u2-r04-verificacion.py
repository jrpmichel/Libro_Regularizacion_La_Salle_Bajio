# ID: ALG-U2-R04
# Libro: Alg. U2 · problema resuelto ALG-U2-04 (verificacion)
# Repositorio: alg/u2_matrices/fragmentos/alg-u2-r04-verificacion.py
import numpy as np

A = np.array([[1.2, 0.4], [0.0, 0.8]])
B = np.array([[1.0, 2.0], [0.5, 1.0]])
C = np.array([[-1.0, 0.0], [0.0, 1.0]])
for nombre, T in (("A", A), ("B", B), ("C", C)):
    print(nombre, "det =", round(np.linalg.det(T), 12),
          " imagenes de i y j:", T[:, 0], T[:, 1])
area = np.linalg.det(A) * 200 * 100
print("area de 200 x 100 px bajo A:", round(area), "px2")     # 19200
# columnas proporcionales: todo el plano cae en una recta
print("columnas de B proporcionales:", np.allclose(B[:, 1], 2 * B[:, 0]))
