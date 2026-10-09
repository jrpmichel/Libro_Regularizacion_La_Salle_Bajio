# ID: DIF-U4-L1
# Libro: Dif. U4 · Laboratorio del Error 4.1 (prueba con código)
# Repositorio: dif/u4_derivada/fragmentos/dif-u4-l1-prueba.py
import sympy as sp

f = lambda x: x**3
for h in [0.1, 0.01, 0.001]:                  # h fija: una secante
    print(h, round((f(1 + h) - f(1)) / h, 6))  # 3.31  3.0301  3.003001

h = sp.symbols("h", real=True)
cociente = sp.expand((f(1 + h) - f(1)) / h)
print(cociente)                                # h**2 + 3*h + 3
print(sp.limit(cociente, h, 0))                # 3: pendiente tangente
