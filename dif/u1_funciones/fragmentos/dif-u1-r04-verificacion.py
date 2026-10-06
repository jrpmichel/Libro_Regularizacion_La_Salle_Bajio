# ID: DIF-U1-R04
# Libro: Dif. U1 · problema resuelto DIF-U1-04 (verificación)
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-r04-verificacion.py
import sympy as sp

q = sp.symbols("q", real=True)
U = q*(500 - 2*q) - (8000 + 100*q)
print(sp.expand(U))                                 # -2*q**2 + 400*q - 8000
print([round(float(r), 2) for r in sp.solve(U, q)]) # [22.54, 177.46]
print(U.subs(q, 100))                               # 12000
print([U.subs(q, k) for k in (22, 23, 177, 178)])   # [-168, 142, 142, -168]
