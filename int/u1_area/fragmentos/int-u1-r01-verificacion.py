# ID: INT-U1-R01
# Libro: Int. U1 · problema resuelto INT-U1-01 (verificación)
# Repositorio: int/u1_area/fragmentos/int-u1-r01-verificacion.py
import numpy as np

def area_numerica(f, a, b, n=100000):    # rectángulos muy delgados
    x = np.linspace(a, b, n + 1)
    return np.sum(f((x[:-1] + x[1:]) / 2)) * (b - a) / n

casos = [("a) rectángulo", lambda x: 4 + 0 * x, 0, 3, 3 * 4),
         ("b) triángulo", lambda x: x / 2, 0, 4, 4 * 2 / 2),
         ("c) cuarto de círculo", lambda x: np.sqrt(4 - x**2), 0, 2,
          np.pi * 2**2 / 4)]
for nombre, f, a, b, formula in casos:
    num = area_numerica(f, a, b)
    print(f"{nombre}: fórmula {formula:.6f}, numérica {num:.6f}")
    assert abs(num - formula) < 1e-6       # la verificación falla si no coinciden
