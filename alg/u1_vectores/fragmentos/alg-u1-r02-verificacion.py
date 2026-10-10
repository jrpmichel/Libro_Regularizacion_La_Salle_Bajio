# ID: ALG-U1-R02
# Libro: Alg. U1 · problema resuelto ALG-U1-02 (verificacion)
# Repositorio: alg/u1_vectores/fragmentos/alg-u1-r02-verificacion.py
import numpy as np

w = np.array([5.0, 0.0])              # viento hacia el este, m/s
a1 = np.array([0.0, 12.0])            # (a) el dron apunta al norte
g1 = a1 + w
print("(a) suelo:", g1, " rapidez:", np.linalg.norm(g1))     # 13.0
desvio = np.degrees(np.arctan2(g1[0], g1[1]))  # desde el norte: (este, norte)
print("    desvio al este:", round(desvio, 2), "grados")      # 22.62

a2 = np.array([-5.0, np.sqrt(12**2 - 5**2)])   # (b) apunta contra el viento
g2 = a2 + w
print("(b) apunta:", np.round(a2, 3), " |a2| =", np.linalg.norm(a2))
print("    suelo:", np.round(g2, 3))
print("    rapidez:", round(np.linalg.norm(g2), 3))               # 10.909
print("    al oeste del norte:", round(np.degrees(np.arcsin(5 / 12)), 2))
print("unitario de (a):", np.round(g1 / np.linalg.norm(g1), 4))  # 5/13, 12/13
