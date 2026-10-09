# ID: INT-U5-S1
# Libro: Int. U5 · subtema 5.1: la sustitución u = x^2 conserva el área
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-s1-sustitucion.py
import numpy as np
import sympy as sp


def punto_medio(f, a, b, n=1000):
    x = np.linspace(a, b, n + 1)
    return np.sum(f((x[:-1] + x[1:]) / 2)) * (b - a) / n


b = np.sqrt(np.pi / 2)
en_x = punto_medio(lambda x: 2 * x * np.cos(x**2), 0, b)
en_u = punto_medio(np.cos, 0, np.pi / 2)        # límites nuevos: 0 y pi/2
print("en x:", round(en_x, 5), "  en u:", round(en_u, 5))   # 1.0 y 1.0

x = sp.symbols("x")
print("(sen x^2)' =", sp.diff(sp.sin(x**2), x))  # 2*x*cos(x**2)
