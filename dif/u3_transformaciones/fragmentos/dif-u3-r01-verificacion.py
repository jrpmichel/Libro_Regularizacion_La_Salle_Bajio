# ID: DIF-U3-R01
# Libro: Dif. U3 · problema resuelto DIF-U3-01 (verificación)
# Repositorio: dif/u3_transformaciones/fragmentos/dif-u3-r01-verificacion.py
import numpy as np

f = lambda t: np.where(t >= 0, 50 * (1 - np.exp(-np.maximum(t, 0) / 5)), 0.0)
g = lambda t: f(t - 2) + 20             # retardo de 2 s, inicio en 20 °C

print("t = 0, 1, 2 s:", g(np.array([0.0, 1.0, 2.0])))       # 20 20 20
print("g(7) =", round(float(g(7.0)), 2), "| f(5) + 20 =",
      round(float(f(5.0)) + 20, 2))                          # 51.61
print("valor final:", round(float(g(200.0)), 4))            # 70.0
