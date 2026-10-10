# ID: ALG-U2-R02
# Libro: Alg. U2 · problema resuelto ALG-U2-02 (verificacion)
# Repositorio: alg/u2_matrices/fragmentos/alg-u2-r02-verificacion.py
import sympy as sp

A = sp.Matrix([[2, -1, 3], [0, 4, 1], [5, 2, -2]])
(a, b, c), (d, e, f), (g, h, i) = A.tolist()
sarrus = a*e*i + b*f*g + c*d*h - c*e*g - a*f*h - b*d*i
col1 = 2 * A.minor(0, 0) - 0 * A.minor(1, 0) + 5 * A.minor(2, 0)
print("Sarrus:", sarrus)                         # -85
print("menores M11, M31:", A.minor(0, 0), A.minor(2, 0))   # -10 y -13
print("cofactores por la columna 1:", col1)      # -85
print("sympy:", A.det())
