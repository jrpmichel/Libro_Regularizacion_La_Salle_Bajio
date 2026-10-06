# ID: DIF-U1-S1
# Libro: Dif. U1 · subtema 1.1: tabla, diferencias y gráfica
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-s1-tabla-a-curva.py
import numpy as np
import matplotlib.pyplot as plt

f = lambda x: x**2 + 1            # regla de asignación
x = np.arange(0, 5)               # entradas: 0, 1, 2, 3, 4
y = f(x)                          # salidas
print("x =", x)
print("y =", y)
print("primera diferencia:", np.diff(y))      # [1 3 5 7]
print("segunda diferencia:", np.diff(y, 2))   # [2 2 2] -> cuadrática
plt.plot(x, y, "o"); plt.grid(True)
plt.xlabel("x"); plt.ylabel("y"); plt.show()
