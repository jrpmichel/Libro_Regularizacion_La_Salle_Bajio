# ID: ALG-U4-R04
# Libro: Alg. U4 · problema resuelto ALG-U4-04 (verificacion)
# Repositorio: alg/u4_espacios/fragmentos/alg-u4-r04-verificacion.py
from itertools import combinations
import sympy as sp

d = [sp.Matrix(v) for v in ([1, 1, 0], [1, 0, 1], [0, 1, 1], [1, 2, 1])]
D = sp.Matrix.hstack(*d)
print("rango de los cuatro:", D.rank())                  # 3: cubren R^3
print("relacion:", list(D.nullspace()[0]))               # d4 = d1 + d3
for trio in combinations(range(4), 3):
    falla = (set(range(4)) - set(trio)).pop() + 1
    det = sp.Matrix.hstack(*[d[i] for i in trio]).det()
    print(f"falla el {falla}: det = {det}", "(base)" if det != 0 else "(NO es base)")
n = d[0].cross(d[2])
print("plano sin el 2, normal:", list(n))                # (1, -1, 1)
print("n . d1, n . d3, n . d4:", [n.dot(d[i]) for i in (0, 2, 3)])   # todos 0
