# ID: ALG-U3-L2
# Libro: Alg. U3 · Laboratorio del Error 3.2 (prueba con codigo)
# Repositorio: alg/u3_sistemas/fragmentos/alg-u3-l2-prueba.py
import numpy as np
import sympy as sp

A = [[2, 1, 3], [1, 2, 1], [3, 3, 4]]      # corte, ensamble, inspeccion
b = [50, 40, 90]                           # horas disponibles
try:
    x = np.linalg.solve(np.array(A, float), b)
    print("numpy devolvio", x, "sin aviso: revisa el residuo")
except np.linalg.LinAlgError as e:
    print("numpy:", e)                     # Singular matrix
M, v = sp.Matrix(A), sp.Matrix(b)
print("rangos:", M.rank(), M.row_join(v).rank())    # 2 2
x1, x2, x3 = sp.symbols("x1:4")
print(sp.linsolve((M, v), x1, x2, x3))              # x3 libre
for t in range(0, 13, 3):                           # lotes enteros
    s = sp.Matrix([20 - sp.Rational(5, 3) * t, 10 + sp.Rational(t, 3), t])
    print(list(s), "usa", list(M * s), "h")         # 50, 40, 90
