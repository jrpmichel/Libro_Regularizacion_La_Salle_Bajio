# ID: DIF-U3-L3
# Libro: Dif. U3 · Laboratorio del Error 3.3 (prueba con código)
# Repositorio: dif/u3_transformaciones/fragmentos/dif-u3-l3-prueba.py
import math

for x in (-4, -1, 0):
    print(f"x = {x}: sqrt(-x) = {math.sqrt(-x)}")   # x <= 0, valores >= 0
try:
    math.sqrt(-4.0)                          # sqrt(-x) en x = 4
except ValueError as error:
    print("x = 4:", error)
print([-math.sqrt(x) for x in (1, 4, 9)])   # -sqrt(x): otra función
