# ID: DIF-U1-L4
# Libro: Dif. U1 · Laboratorio del Error 1.4 (prueba con código)
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-l4-prueba.py
import math

L = 0.80                                              # m, longitud del brazo
print(round(L * math.sin(90), 4))                     # 0.7152: math.sin recibe radianes
print(round(L * math.sin(math.radians(90)), 4))       # 0.8: el valor correcto
print(math.isclose(math.sin(math.radians(90)), 1.0),  # True
      math.isclose(math.sin(90), 1.0))                # False: la prueba del caso límite lo delata
