# ID: ALG-U1-R04
# Libro: Alg. U1 · problema resuelto ALG-U1-04 (verificacion)
# Repositorio: alg/u1_vectores/fragmentos/alg-u1-r04-verificacion.py
import numpy as np

r = np.array([1.5, 0.0, 6.0])          # m, de la base al centro del anuncio
F = np.array([0.0, 2400.0, -1500.0])   # N, viento + peso
M = np.cross(r, F)
print("M = r x F =", M, "N.m")                         # [-14400 2250 3600]
print("|M| =", round(np.linalg.norm(M), 1), "N.m")     # 15012.7
print("M . r =", np.dot(M, r), " M . F =", np.dot(M, F))   # 0 y 0
print("F x r =", np.cross(F, r))       # el orden invertido cambia el signo
