# ID: INT-U5-R02
# Libro: Int. U5 · problema resuelto INT-U5-02 (verificación)
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-r02-verificacion.py
import numpy as np
import sympy as sp

m, L, g, T = 2, 0.3, 9.81, 2.0
u = sp.symbols("u")
W = m * g * L * sp.integrate(sp.sin(u), (u, 0, sp.pi / 2))   # límites en u
print("W =", round(float(W), 3), "J")                        # 5.886

n = 200
t = (np.arange(n) + 0.5) * T / n                    # puntos medios
th = np.pi / 2 * (3 * (t / T)**2 - 2 * (t / T)**3)  # theta(t)
dth = np.pi / 2 * (6 * t / T**2 - 6 * t**2 / T**3)  # theta'(t)
p = m * g * L * np.sin(th) * dth                    # potencia, W
print("punto medio con p(t):", round(p.sum() * T / n, 3), "J")   # 5.886
