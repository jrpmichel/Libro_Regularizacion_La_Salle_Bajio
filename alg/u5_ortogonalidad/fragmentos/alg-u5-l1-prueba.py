# ID: ALG-U5-L1
# Libro: Alg. U5 · Laboratorio del Error 5.1 (prueba con codigo)
# Repositorio: alg/u5_ortogonalidad/fragmentos/alg-u5-l1-prueba.py
import numpy as np

g1, g2 = np.array([1, 0, 1]), np.array([0, 1, 1])   # guias de la mesa
v = np.array([2, 2, 1])                              # velocidad mandada (cm/s)
proy = lambda x, u: (x @ u) / (u @ u) * u
suma = proy(v, g1) + proy(v, g2)                     # error: base no ortogonal
print("g1.g2 =", g1 @ g2, "  suma =", suma)
print("|suma| =", round(np.linalg.norm(suma), 3), "> |v| =", np.linalg.norm(v))
print("(v - suma).g1 =", (v - suma) @ g1)            # no es perpendicular
A = np.column_stack([g1, g2])
p = A @ np.linalg.solve(A.T @ A, A.T @ v)           # ecuaciones normales
n = np.cross(g1, g2)
print("correcta p =", p, "  con la normal:", v - (v @ n) / (n @ n) * n)
print("|p| =", round(np.linalg.norm(p), 3))
print("(v - p).g1 =", (v - p) @ g1, "  (v - p).g2 =", (v - p) @ g2)
