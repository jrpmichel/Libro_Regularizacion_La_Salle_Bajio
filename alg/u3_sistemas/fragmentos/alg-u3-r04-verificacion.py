# ID: ALG-U3-R04
# Libro: Alg. U3 · problema resuelto ALG-U3-04 (verificacion)
# Repositorio: alg/u3_sistemas/fragmentos/alg-u3-r04-verificacion.py
import sympy as sp

Ab = sp.Matrix([[1, 0, 0, -1, 50], [-1, 1, 0, 0, -40],
                [0, -1, 1, 0, 30], [0, 0, -1, 1, -40]])
A, b = Ab[:, :4], Ab[:, 4]
print("rangos:", A.rank(), Ab.rank())               # 3 3
x1, x2, x3, x4 = sp.symbols("x1:5")
print(sp.linsolve((A, b), x1, x2, x3, x4))          # x4 libre
t = sp.symbols("t")
x = [50 + t, 10 + t, 40 + t, t]                     # tarimas/h
cond = [xi >= 0 for xi in x] + [xi <= 60 for xi in x]
print("intervalo de t:", sp.reduce_inequalities(cond, t))   # 0..10
b2 = b.copy()
b2[2] = 35                                          # entran 35 en C
print("con 35 en C:", A.rank(), A.row_join(b2).rank())      # 3 4
