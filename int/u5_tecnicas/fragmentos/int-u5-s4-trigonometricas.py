# ID: INT-U5-S4
# Libro: Int. U5 · subtema 5.4: área de sen^2 y del cuarto de círculo
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-s4-trigonometricas.py
import numpy as np


def punto_medio(f, a, b, n=2000):
    x = np.linspace(a, b, n + 1)
    return np.sum(f((x[:-1] + x[1:]) / 2)) * (b - a) / n


s2 = punto_medio(lambda x: np.sin(x)**2, 0, 2 * np.pi)
print("sen^2 en [0, 2pi]:", round(s2, 5), "  pi =", round(np.pi, 5))
a = 1.0
c = punto_medio(lambda x: np.sqrt(a**2 - x**2), 0, a)
# con x = a sen(t) queda a^2 cos^2(t) en [0, pi/2], sin raíz
t = punto_medio(lambda t: a**2 * np.cos(t)**2, 0, np.pi / 2)
print("cuarto de círculo:", round(c, 5), " sustituida:", round(t, 5),
      " pi/4 =", round(np.pi / 4, 5))
