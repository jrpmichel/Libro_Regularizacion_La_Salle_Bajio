# ID: INT-U6-R03
# Libro: Int. U6 · problema resuelto INT-U6-03 (verificación)
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-r03-verificacion.py
import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

def exacta(s, d):
    """Longitud del cable y = 4 d x^2 / s^2 entre x = -s/2 y s/2."""
    f = lambda x: np.sqrt(1 + (8 * d * x / s**2)**2)
    return quad(f, -s / 2, s / 2)[0]

s = 200
for d in (20, 50):
    L, M = exacta(s, d), s * (1 + 8 / 3 * (d / s)**2)
    print(f"d = {d} m: exacta {L:.3f}, manual {M:.3f},"
          f" error {100 * (M - L) / L:.3f} %")
n = brentq(lambda n: (1 + 8 / 3 * n**2) / exacta(1, n) - 1.005, 0.05, 0.4)
print("el error llega a 0.5 % con d/s =", round(n, 3))   # 0.177
