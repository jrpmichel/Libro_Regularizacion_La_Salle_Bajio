# ID: DIF-U5-L3
# Libro: Dif. U5 · Laboratorio del Error 5.3 (prueba con código)
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-l3-prueba.py
import math

f = lambda x: x**x
escena = lambda x: x * x**(x - 1)               # la respuesta: x^x
bien = lambda x: x**x * (math.log(x) + 1)       # derivación logarítmica
h = 1e-6
for a in (1.0, 2.0, 0.5):
    c = (f(a + h) - f(a)) / h
    print(f"x = {a}: cociente {c:.4f} | escena {escena(a):.4f}"
          f" | correcta {bien(a):.4f}")
