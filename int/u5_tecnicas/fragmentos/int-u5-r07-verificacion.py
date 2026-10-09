# ID: INT-U5-R07
# Libro: Int. U5 · problema resuelto INT-U5-07 (verificación)
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-r07-verificacion.py
import numpy as np
from scipy.integrate import quad

r, L = 1.0, 4.0


def A(h):                          # fórmula obtenida con y = r sen(t)
    y = h - r
    return (r**2 * np.arcsin(y / r) + y * np.sqrt(r**2 - y**2)
            + np.pi * r**2 / 2)


num = quad(lambda y: 2 * np.sqrt(r**2 - y**2), -r, 0.5 - r)[0]
print("A(0.5):", round(A(0.5), 4), " cuadratura:", round(num, 4))  # 0.6142
V = L * A(0.5)
print("V =", round(V, 3), "m3,", round(100 * V / (np.pi * r**2 * L), 1), "%")
print("A(r) =", round(A(r), 6), " A(2r) =", round(A(2 * r), 6))  # pi/2, pi
