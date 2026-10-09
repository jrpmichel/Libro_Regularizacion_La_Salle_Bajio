# ID: DIF-U2-R01
# Libro: Dif. U2 · problema resuelto DIF-U2-01 (verificación)
# Repositorio: dif/u2_limite/fragmentos/dif-u2-r01-verificacion.py
import numpy as np
import sympy as sp

factor = lambda i, n=12: (1 - (1 + i)**-n) / i     # no existe en i = 0
for i in (0.01, 0.001, 0.0001):
    print(f"i = {i:<7} {factor(i):.4f}   |   i = {-i:<8} {factor(-i):.4f}")

j = sp.symbols("i")
print("límite:", sp.limit((1 - (1 + j)**-12) / j, j, 0))       # 12
suma = lambda i: sum((1 + i)**-k for k in range(1, 13))      # como suma
print(np.isclose(suma(0.01), factor(0.01)), suma(0))         # True 12.0
