# ID: ALG-U1-R03
# Libro: Alg. U1 · problema resuelto ALG-U1-03 (verificacion)
# Repositorio: alg/u1_vectores/fragmentos/alg-u1-r03-verificacion.py
import numpy as np

F = np.array([500.0, 120.0, 450.0])   # N
d = np.array([8.0, 0.0, 6.0])         # m, de la base a la cima de la rampa
W = np.dot(F, d)
nF = np.linalg.norm(F)
t = np.degrees(np.arccos(np.clip(W / (nF * np.linalg.norm(d)), -1, 1)))
u = d / np.linalg.norm(d)
Fpar = np.dot(F, u)
Fperp = F - Fpar * u
print("W =", W, "J")                                       # 6700
print("|F| =", round(nF, 1), "N   angulo =", round(t, 1), "grados")
print("F paralela =", Fpar, "N")                           # 670
print("F perpendicular =", Fperp, "N")                     # [-36 120 48]
print("|F perpendicular| =", round(np.linalg.norm(Fperp), 1), "N")
print("Fpar * |d| =", Fpar * np.linalg.norm(d))            # otra vez 6700
