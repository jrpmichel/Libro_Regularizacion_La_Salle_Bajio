# ID: ALG-U3-R05
# Libro: Alg. U3 · problema resuelto ALG-U3-05 (verificacion)
# Repositorio: alg/u3_sistemas/fragmentos/alg-u3-r05-verificacion.py
import sympy as sp

x, y, z = sp.symbols("x y z")
A = sp.Matrix([[1, 1, 1], [1, 2, 3], [2, 1, 1]])
print("(a)", sp.linsolve((A, sp.Matrix([6, 14, 7])), x, y, z))
B = sp.Matrix([[1, 1, 1], [1, 2, 3], [3, 4, 5]])
for c in (26, 30):                                  # (b) y (c)
    bc = sp.Matrix([6, 14, c])
    print(c, "rangos:", B.rank(), B.row_join(bc).rank(),
          sp.linsolve((B, bc), x, y, z))
n = [sp.Matrix(f) for f in B.tolist()]              # normales
pares = ((0, 1), (0, 2), (1, 2))
print("normales paralelas:", [n[i].cross(n[j]) == sp.zeros(3, 1)
                              for i, j in pares])   # todas False
