# ID: ALG-U1-S5
# Libro: Alg. U1 · subtema 1.5: plano por tres puntos, distancia e interseccion con una recta
# Repositorio: alg/u1_vectores/fragmentos/alg-u1-s5-recta-plano.py
import numpy as np

A = np.array([1.0, 2.0, 0.0])
B = np.array([3.0, 0.0, 1.0])
C = np.array([0.0, 1.0, 2.0])
n = np.cross(B - A, C - A)            # normal al plano
d = np.dot(n, A)                      # ecuacion: n . (x, y, z) = d
print("normal =", n, " d =", d)       # [-3. -5. -4.]  -13: 3x + 5y + 4z = 13

P = np.array([2.0, 3.0, 4.0])
dist = abs(np.dot(n, P) - d) / np.linalg.norm(n)   # dividir entre |n|
print("distancia de P al plano =", round(dist, 4))

P0 = np.array([0.0, 0.0, 0.0])        # recta: P0 + t u
u = np.array([1.0, 1.0, 1.0])
den = np.dot(n, u)
if abs(den) < 1e-12:
    print("recta paralela al plano")
else:
    t = (d - np.dot(n, P0)) / den
    print("t =", round(t, 4), " punto =", np.round(P0 + t * u, 4))
