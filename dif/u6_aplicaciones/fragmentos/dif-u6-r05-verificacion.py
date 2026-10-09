# ID: DIF-U6-R05
# Libro: Dif. U6 · problema resuelto DIF-U6-05 (verificación)
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-r05-verificacion.py
import numpy as np

L = lambda w: w**4 - 4 * w**2 + w
dL = lambda w: 4 * w**3 - 8 * w + 1
d2L = lambda w: 12 * w**2 - 8


def newton(w, pasos=30):                     # resuelve L'(w) = 0
    for _ in range(pasos):
        w = w - dL(w) / d2L(w)
    return w


for w0 in (-1.5, 0.0, 1.5):
    c = newton(w0)
    print(f"w = {c:.4f}  L = {L(c):.4f}  L'' = {d2L(c):.3f}")
print("L'(0.5) =", dL(0.5), "| L'(0) =", dL(0))   # signo: hacia dónde baja
print("np.roots:", np.sort(np.roots([4, 0, -8, 1])).round(4))   # control
