# ID: INT-U1-S1
# Libro: Int. U1 · subtema 1.1: área bajo una poligonal con trapecios
# Repositorio: int/u1_area/fragmentos/int-u1-s1-poligonal.py
import numpy as np
import matplotlib.pyplot as plt

x = np.array([0, 1, 3, 4])          # vértices de la poligonal (figura 1.1c)
y = np.array([1, 3, 3, 0])
anchos = np.diff(x)                 # base de cada trapecio
alturas = (y[:-1] + y[1:]) / 2      # promedio de los lados paralelos
areas = anchos * alturas
for k, a in enumerate(areas, 1):
    print(f"trapecio {k}: {a}")
print("área total:", areas.sum())   # 2 + 6 + 1.5 = 9.5

plt.fill_between(x, 0, y, color="0.85")
plt.plot(x, y, "o-")
plt.xlabel("x"); plt.ylabel("y")
plt.show()
