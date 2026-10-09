# ID: INT-U5-S2
# Libro: Int. U5 · subtema 5.2: dos áreas que llenan el rectángulo e x 1
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-s2-partes.py
import numpy as np


def punto_medio(f, a, b, n=1000):
    x = np.linspace(a, b, n + 1)
    return np.sum(f((x[:-1] + x[1:]) / 2)) * (b - a) / n


gris = punto_medio(np.log, 1, np.e)      # bajo y = ln x, de 1 a e
rayada = punto_medio(np.exp, 0, 1)       # a la izquierda: x = e^y, y de 0 a 1
print("gris:", round(gris, 5), " rayada:", round(rayada, 5))   # 1.0 y 1.71828
print("suma:", round(gris + rayada, 5), " e =", round(np.e, 5))  # 2.71828
