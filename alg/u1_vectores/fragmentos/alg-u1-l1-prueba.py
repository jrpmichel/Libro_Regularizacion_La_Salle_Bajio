# ID: ALG-U1-L1
# Libro: Alg. U1 · Laboratorio del Error 1.1 (prueba con codigo)
# Repositorio: alg/u1_vectores/fragmentos/alg-u1-l1-prueba.py
import numpy as np

P = np.array([0.0, 0.0, 12.0])         # punta del mastil, m
Q = np.array([4.0, 3.0, 0.0])          # ancla en el suelo, m
T = 1300.0                             # N
mal = T * (P - Q) / np.linalg.norm(P - Q)    # del ancla a la punta
bien = T * (Q - P) / np.linalg.norm(Q - P)   # de la punta hacia el ancla
print("mal: ", mal, " |.| =", np.linalg.norm(mal))
print("bien:", bien, " |.| =", np.linalg.norm(bien))
print("el cable jala hacia abajo:", bien[2] < 0)          # True
print("compresion extra en el mastil:", -bien[2], "N")      # 1200
