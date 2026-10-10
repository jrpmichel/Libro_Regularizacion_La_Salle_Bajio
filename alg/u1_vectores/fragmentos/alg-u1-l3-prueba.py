# ID: ALG-U1-L3
# Libro: Alg. U1 · Laboratorio del Error 1.3 (prueba con codigo)
# Repositorio: alg/u1_vectores/fragmentos/alg-u1-l3-prueba.py
import numpy as np

n, d0 = np.array([1.0, 2.0, 2.0]), 12.0      # talud: x + 2y + 2z = 12
P = np.array([4.0, 5.0, 6.0])                # dron, m
s = np.dot(n, P) - d0
print("respuesta de la IA:", abs(s), "m")                          # 14
dist = abs(s) / np.linalg.norm(n)
print("distancia:", round(dist, 4), "m")                           # 4.6667
H = P - (s / np.dot(n, n)) * n               # pie de la perpendicular
print("H =", np.round(H, 4), " en el plano:", np.isclose(np.dot(n, H), d0))
print("|PH| =", round(np.linalg.norm(P - H), 4), "m  respeta 5 m:", dist >= 5)
