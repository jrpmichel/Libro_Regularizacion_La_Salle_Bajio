# ID: ALG-U6-S1
# Libro: Alg. U6 · subtema 6.1: tensiones de tres cables, momento y trabajo
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-s1-equilibrio.py
import sympy as sp

anclajes = [sp.Matrix(a) for a in ([2, 2, 1], [-2, 1, 2], [1, -2, 2])]  # m
U = sp.Matrix.hstack(*[a / a.norm() for a in anclajes])  # columnas u1 u2 u3
W = sp.Matrix([0, 0, -90])                         # peso de la lampara (N)
T = U.solve(-W)                                    # U T = -W
print("tensiones (N):", list(T))                   # 30, 60, 60
print("U^T U = I:", U.T * U == sp.eye(3))          # columnas ortonormales
r = sp.Matrix([sp.Rational(3, 10), sp.Rational(4, 10), 0])  # m
F = sp.Matrix([0, 0, -200])                                # N
M = r.cross(F)                                     # r x F, en ese orden
print("M =", list(M), "N m;  |M| =", M.norm())
print("F x r =", list(F.cross(r)), "(el signo cambia)")
F2, d = sp.Matrix([40, 30, 0]), sp.Matrix([5, 0, 0])
print("trabajo W = F . d =", F2.dot(d), "J")
