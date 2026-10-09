# ID: INT-U2-S4
# Libro: Int. U2 · subtema 2.4: propiedades de la integral comprobadas con sumas
# Repositorio: int/u2_riemann/fragmentos/int-u2-s4-propiedades.py
import numpy as np

def I(f, a, b, n=100000):            # suma de punto medio; acepta a > b
    x = np.linspace(a, b, n + 1)
    return np.sum(f((x[:-1] + x[1:]) / 2)) * (b - a) / n

f = lambda x: x**2
g = lambda x: np.exp(-x**2)
print("linealidad:", round(I(lambda x: 3 * f(x) - 2 * x + 5, 0, 1), 6))   # 5
print("aditividad:", round(I(f, 1, 2), 6), round(I(f, 0, 2) - I(f, 0, 1), 6))  # 7/3
print("límites al revés:", round(I(f, 2, 1), 6))                    # -7/3
print("cotas:", round(np.exp(-1), 4), "<=", round(I(g, 0, 1), 4), "<= 1")
print("impar:", round(I(lambda x: x**3, -2, 2), 6), "  par:", round(I(f, -2, 2), 6))
