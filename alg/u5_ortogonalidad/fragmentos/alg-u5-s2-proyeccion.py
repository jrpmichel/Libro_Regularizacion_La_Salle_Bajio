# ID: ALG-U5-S2
# Libro: Alg. U5 · subtema 5.2: proyeccion sobre un plano con las ecuaciones normales
# Repositorio: alg/u5_ortogonalidad/fragmentos/alg-u5-s2-proyeccion.py
import sympy as sp

A = sp.Matrix([[1, 1], [0, 1], [1, 0]])   # columnas a1 = (1,0,1), a2 = (1,1,0)
b = sp.Matrix([1, 2, 3])
c = (A.T * A).solve(A.T * b)               # ecuaciones normales A^T A c = A^T b
p = A * c                                  # punto del plano mas cercano a b
r = b - p
print("c =", list(c), " p =", list(p), " r =", list(r))
print("A^T r =", list(A.T * r))            # [0, 0]: r perpendicular al plano
print("distancia de b al plano:", r.norm(), "=", round(float(r.norm()), 4))
w = A * sp.Matrix([1, 1])                  # otro punto del plano
print("|b - w| =", round(float((b - w).norm()), 4), "> |r|")
P = A * (A.T * A).inv() * A.T              # matriz de la proyeccion
print("P^2 == P:", P * P == P, "  P^T == P:", P.T == P)
