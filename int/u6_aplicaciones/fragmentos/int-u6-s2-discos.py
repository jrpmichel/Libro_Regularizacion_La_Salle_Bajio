# ID: INT-U6-S2
# Libro: Int. U6 · subtema 6.2: volumen de un cono sumando discos
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-s2-discos.py
import numpy as np

r, h = 3.0, 4.0                         # radio y altura del cono (m)
radio = lambda x: r / h * x             # radio del disco en x
for n in (4, 16, 64, 256):
    dx = h / n
    xm = (np.arange(n) + 0.5) * dx      # puntos medios
    V = np.sum(np.pi * radio(xm)**2 * dx)
    print(f"n = {n:3d}: V = {V:.4f} m^3")
print("exacto, pi r^2 h / 3 =", round(np.pi * r**2 * h / 3, 4))  # 37.6991
