# ID: ALG-U3-L1
# Libro: Alg. U3 · Laboratorio del Error 3.1 (prueba con codigo)
# Repositorio: alg/u3_sistemas/fragmentos/alg-u3-l1-prueba.py
from fractions import Fraction as F
import numpy as np

A = np.array([[1, 1, 1], [2, 3, 1], [1, 2, 3]])     # minutos por bomba
b = np.array([12, 25, 25])                          # litros
for nombre, b3 in (("con error", 13), ("correcto", 12)):
    z = F(b3, 3)                                    # renglon 3: 3z = b3
    y = 1 + z                                       # renglon 2: y - z = 1
    x = 12 - y - z                                  # renglon 1
    s = np.array([x, y, z])
    residuo = b - A @ s
    print(nombre, [str(v) for v in s],
          "residuo:", [str(r) for r in residuo])    # solo falla el 3
