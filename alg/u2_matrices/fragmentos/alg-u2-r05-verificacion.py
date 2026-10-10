# ID: ALG-U2-R05
# Libro: Alg. U2 · problema resuelto ALG-U2-05 (verificacion)
# Repositorio: alg/u2_matrices/fragmentos/alg-u2-r05-verificacion.py
import sympy as sp

K = sp.Matrix([[4, -2, 0], [-2, 4, -2], [0, -2, 2]])   # kN/mm
L, U, _ = K.LUdecomposition()                          # sin intercambios
print("L =", L.tolist())
print("U =", U.tolist())
print("det K =", K.det(), "=", U[0, 0] * U[1, 1] * U[2, 2])
for F in ([10, 20, 30], [0, 0, 10], [5, 5, 5]):        # kN
    y = L.solve(sp.Matrix(F))                          # hacia adelante
    x = U.solve(y)                                     # hacia atras
    print("F =", F, " y =", list(y), " x =", list(x), "mm")
