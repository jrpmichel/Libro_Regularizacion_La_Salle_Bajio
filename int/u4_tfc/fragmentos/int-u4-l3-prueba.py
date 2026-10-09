# ID: INT-U4-L3
# Libro: Int. U4 · Laboratorio del Error 4.3 (prueba con código)
# Repositorio: int/u4_tfc/fragmentos/int-u4-l3-prueba.py
import numpy as np

grados = np.arange(0, 91, 5.0)             # columna de la hoja
tau = 12 * np.cos(np.radians(grados))      # N m
def trapecio(y, xs):
    return float(np.sum((y[:-1] + y[1:]) / 2 * np.diff(xs)))
print("ancho en grados:", round(trapecio(tau, grados), 1), "J")   # 687
print("ancho en radianes:", round(trapecio(tau, np.radians(grados)), 2))
print("cota 12 * pi/2 =", round(12 * np.pi / 2, 1), "J")          # 18.8
