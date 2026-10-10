# ID: ALG-U3-S4
# Libro: Alg. U3 · subtema 3.4: clasificacion con rangos y solucion general
# Repositorio: alg/u3_sistemas/fragmentos/alg-u3-s4-clasificacion.py
import sympy as sp

def clasifica(Ab):
    Ab = sp.Matrix(Ab)
    A, b, n = Ab[:, :-1], Ab[:, -1], Ab.cols - 1
    r, rs = A.rank(), Ab.rank()
    if r < rs:
        return f"incompatible (r = {r}, r* = {rs})"
    xs = sp.symbols(f"x1:{n + 1}")
    sol = sp.linsolve((A, b), *xs)
    if r == n:
        return f"determinado: {sol}"
    return f"indeterminado, {n - r} parametro(s): {sol}"

print(clasifica([[1, 1, 1, 6], [2, 2, 3, 15], [1, 3, 2, 13]]))
print(clasifica([[1, 2, -1, 3], [2, 4, 1, 9]]))      # x2 libre
print(clasifica([[1, 1, 1, 1], [2, 2, 2, 3], [1, -1, 0, 0]]))
H = sp.Matrix([[3, 0, -1, 0], [8, 0, 0, -2], [0, 2, -2, -1]])
v = H.nullspace()[0]                     # propano: a, b, c, d
print("homogeneo:", list(v / v[0]))      # [1, 5, 3, 4]
