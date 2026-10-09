# ID: INT-U6-S5
# Libro: Int. U6 · subtema 6.5: masa, centro de masa e inercia por trozos
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-s5-masa.py
import numpy as np

L, n = 0.5, 1000
dx = L / n
x = (np.arange(n) + 0.5) * dx           # puntos medios (m)
dm = (4 - 4 * x) * dx                   # masa de cada trozo (kg)
m = dm.sum()
xc = (x * dm).sum() / m
I0 = (x**2 * dm).sum()
print(f"m = {m:.4f} kg, x_c = {xc:.4f} m, I_O = {I0:.5f} kg m^2")
print("exactos: 1.5,", round(2 / 9, 4), "y", round(5 / 48, 5))
print("barra uniforme de 1.5 kg, m L^2/3 =", round(1.5 * L**2 / 3, 5))
