# ID: INT-U6-S3
# Libro: Int. U6 · subtema 6.3: poligonales sobre y = x^(3/2)
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-s3-poligonal.py
import numpy as np

f = lambda x: x**1.5
exacta = 8 / 27 * (10**1.5 - 1)                      # 9.07342
for n in (1, 2, 4, 8, 64):
    x = np.linspace(0, 4, n + 1)
    L = np.sum(np.hypot(np.diff(x), np.diff(f(x))))  # suma de segmentos
    print(f"n = {n:2d}: poligonal {L:.5f}, falta {exacta - L:.5f}")
print("integral:", round(exacta, 5))
