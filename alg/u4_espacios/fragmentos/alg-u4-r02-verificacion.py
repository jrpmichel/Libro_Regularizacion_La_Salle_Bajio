# ID: ALG-U4-R02
# Libro: Alg. U4 · problema resuelto ALG-U4-02 (verificacion)
# Repositorio: alg/u4_espacios/fragmentos/alg-u4-r02-verificacion.py
import numpy as np
import sympy as sp

estado = lambda t: 10 * np.cos(np.deg2rad([t, t - 120, t + 120]))   # (ia, ib, ic)
for t in (0, 60, 120):
    print("theta =", t, " i =", np.round(estado(t), 6))
V = sp.Matrix([[10, 5, -5], [-5, 5, 10], [-5, -10, -5]])   # columnas i1, i2, i3
print("rango:", V.rank())                                  # 2: dependientes
print("relacion:", list(V.nullspace()[0]))                 # (1, -1, 1)
n = V[:, 0].cross(V[:, 1])
print("normal:", list(n), "-> plano ia + ib + ic = 0")     # (75, 75, 75)
print("suma en 500 instantes:", np.abs(estado(np.linspace(0, 360, 500)).sum(axis=0)).max())
print("ic con 7.3 y -9.6 A:", round(-(7.3 - 9.6), 2), "A")
