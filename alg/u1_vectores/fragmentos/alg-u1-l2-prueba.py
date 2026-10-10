# ID: ALG-U1-L2
# Libro: Alg. U1 · Laboratorio del Error 1.2 (prueba con codigo)
# Repositorio: alg/u1_vectores/fragmentos/alg-u1-l2-prueba.py
import numpy as np

a = np.array([0.30, 0.30, 0.00])       # m, brazo
b = np.array([-0.25, 0.00, 0.25])      # m, antebrazo
na, nb = np.linalg.norm(a), np.linalg.norm(b)
s = np.linalg.norm(np.cross(a, b)) / (na * nb)     # seno
c = np.dot(a, b) / (na * nb)                        # coseno
print("con arcsin:", round(np.degrees(np.arcsin(s)), 2))    # 60: mal
print("con arccos:", round(np.degrees(np.arccos(c)), 2))    # 120
print("con arctan2:", round(np.degrees(np.arctan2(s, c)), 2))   # 120
print("a . b =", round(np.dot(a, b), 4), "-> negativo: angulo obtuso")
