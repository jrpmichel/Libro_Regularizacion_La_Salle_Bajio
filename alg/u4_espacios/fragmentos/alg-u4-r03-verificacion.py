# ID: ALG-U4-R03
# Libro: Alg. U4 · problema resuelto ALG-U4-03 (verificacion)
# Repositorio: alg/u4_espacios/fragmentos/alg-u4-r03-verificacion.py
import numpy as np
import sympy as sp

h = sp.Rational(1, 2)
P = sp.Matrix([[h, h], [h, -h]])          # columnas: mA y mB, mm por mm de banda
print("det P =", P.det(), " P^(-1) =", P.inv().tolist())
ab = P.inv() * sp.Matrix([30, 10])
print("bandas para (30, 10):", list(ab), " regreso:", list(P * ab))
for lado in ([20, 0], [0, 20], [-20, 0], [0, -20]):        # cuadrado de 20 mm
    print("lado", lado, "-> bandas", list(P.inv() * sp.Matrix(lado)))
Pi = np.array(P.inv(), dtype=float)
for nombre, d in (("eje x", [1, 0]), ("diagonal", [1, 1] / np.sqrt(2))):
    vmax = 200 / np.max(np.abs(Pi @ d))   # ninguna banda pasa de 200 mm/s
    print(nombre, "rapidez maxima:", round(vmax, 1), "mm/s")
