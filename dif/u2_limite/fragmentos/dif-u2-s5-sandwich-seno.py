# ID: DIF-U2-S5
# Libro: Dif. U2 · subtema 2.5: sándwich y sen x / x
# Repositorio: dif/u2_limite/fragmentos/dif-u2-s5-sandwich-seno.py
import numpy as np

for x in (0.5, 0.1, 0.01, 0.001):                  # x en radianes
    print(f"x = {x:<6} cos x = {np.cos(x):.8f}  "
          f"sen x / x = {np.sin(x)/x:.8f}")

g = 0.01                                           # x en grados
print("en grados:", np.sin(np.radians(g)) / g, "aprox. pi/180 =", np.pi / 180)

x = np.array([0.1, 0.01, 0.001])
print("x sen(1/x):", np.round(x * np.sin(1 / x), 6), "  cota |x|:", x)
