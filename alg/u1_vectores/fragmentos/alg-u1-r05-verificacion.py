# ID: ALG-U1-R05
# Libro: Alg. U1 · problema resuelto ALG-U1-05 (verificacion)
# Repositorio: alg/u1_vectores/fragmentos/alg-u1-r05-verificacion.py
import numpy as np

n, d0 = np.array([2.0, 1.0, 2.0]), 12.0          # muro: 2x + y + 2z = 12
C = np.array([0.0, 0.0, 1.0])                    # camara
for u in ([1, 1, 1], [1, 0, -1], [-1, -1, 0]):
    u = np.array(u, dtype=float)
    num, den = d0 - np.dot(n, C), np.dot(n, u)
    if abs(den) < 1e-12:
        print(u, "paralelo; num =", num)
    else:
        t = num / den
        print(u, "t =", round(t, 4), " punto =", np.round(C + t * u, 4),
              "visible" if t > 0 else "detras de la camara")
D = abs(np.dot(n, C) - d0) / np.linalg.norm(n)
print("distancia camara-muro =", round(D, 4))      # 3.3333
