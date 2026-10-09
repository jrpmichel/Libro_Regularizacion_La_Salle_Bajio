# ID: INT-U2-S1
# Libro: Int. U2 · subtema 2.1: sumas izquierda, derecha y de punto medio
# Repositorio: int/u2_riemann/fragmentos/int-u2-s1-izq-der-med.py
import numpy as np

def sumas(f, a, b, n):
    x = np.linspace(a, b, n + 1)
    dx = (b - a) / n
    L = np.sum(f(x[:-1])) * dx                 # extremo izquierdo
    R = np.sum(f(x[1:])) * dx                  # extremo derecho
    M = np.sum(f((x[:-1] + x[1:]) / 2)) * dx   # punto medio
    return L, R, M

f = lambda x: x**2                             # valor exacto: 1/3
print("   n        L_n        R_n        M_n")
for n in (4, 10, 100):
    L, R, M = sumas(f, 0, 1, n)
    print(f"{n:4d}  {L:.6f}  {R:.6f}  {M:.6f}")
