# ID: ALG-U2-R01
# Libro: Alg. U2 · problema resuelto ALG-U2-01 (verificacion)
# Repositorio: alg/u2_matrices/fragmentos/alg-u2-r01-verificacion.py
import numpy as np

M = np.array([[0.8, 2.5, 4.0],      # madera, m2 por pieza
              [12, 20, 32],         # tornillos por pieza
              [0.2, 0.6, 0.9]])     # barniz, L por pieza
P = np.array([[40, 25],             # sillas en las semanas 1 y 2
              [10, 15],             # mesas
              [6, 8]])              # libreros
c = np.array([180, 0.5, 120])       # pesos por m2, por tornillo y por L
R = M @ P
print("materiales por semana:\n", R)
print("costo por semana:", c @ R)               # [17344. 19082.]
print("costo por pieza:", c @ M)                # [174. 532. 844.]
print("(cM)P == c(MP):", np.allclose((c @ M) @ P, c @ (M @ P)))
