# ID: INT-U6-R06
# Libro: Int. U6 · problema resuelto INT-U6-06 (verificación)
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-r06-verificacion.py
import math

# distancia real: área bajo la velocidad (trapecio de 2 + 16 + 2 s)
print("distancia real:", 0.5 * 2 * 1 + 16 * 1 + 0.5 * 2 * 1, "m")   # 18
b, tol = 0.05, 0.5                      # sesgo (m/s^2) y tolerancia (m)
for t in (20, 60):
    print(f"t = {t} s: error {0.5 * b * t**2:.0f} m,"
          f" velocidad falsa {b * t:.0f} m/s")
for bb in (0.05, 0.005):
    print(f"b = {bb}: error menor que {tol} m hasta"
          f" t = {math.sqrt(2 * tol / bb):.2f} s")       # 4.47; 14.14
