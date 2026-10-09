# ID: INT-U5-L1
# Libro: Int. U5 · Laboratorio del Error 5.1 (prueba con código)
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-l1-prueba.py
from math import log
from scipy.integrate import quad

mal = 0.5 * log(3)            # límites de x usados con u
bien = 0.5 * log(5)           # u va de 2 a 10
num = quad(lambda x: x / (x**2 + 1), 1, 3)[0]
print("incorrecto:", round(mal, 3), " correcto:", round(bien, 3),
      " cuadratura:", round(num, 3))         # 0.549, 0.805, 0.805
print("cotas: entre", 2 * 0.3, "y", 2 * 0.5)  # 0.6 y 1.0
