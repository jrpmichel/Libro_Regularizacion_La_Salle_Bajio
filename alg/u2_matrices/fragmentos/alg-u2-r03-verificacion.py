# ID: ALG-U2-R03
# Libro: Alg. U2 · problema resuelto ALG-U2-03 (verificacion)
# Repositorio: alg/u2_matrices/fragmentos/alg-u2-r03-verificacion.py
import sympy as sp

A = sp.Matrix([[2, 1, 0], [0, 2, 1], [1, 0, 2]])   # mV/N
v = sp.Matrix([11, 12, 7])                          # mV medidos
print("det A =", A.det())                           # 9
print("adj A =", A.adjugate())
Ainv = A.adjugate() / A.det()
print("A adj A = 9 I:", A * A.adjugate() == 9 * sp.eye(3))
f = Ainv * v
print("f = A^-1 v =", list(f), "N")                 # [3, 5, 2]
print("sin corregir, v/2 =", [float(x) / 2 for x in v], "N")
