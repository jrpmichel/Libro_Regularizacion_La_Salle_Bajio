# ID: ALG-U6-S5
# Libro: Alg. U6 · subtema 6.5: recta por el origen y modelo de primer orden
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-s5-ajuste.py
import numpy as np

I = np.array([0.1, 0.2, 0.3, 0.4])            # A
V = np.array([1.9, 4.1, 6.0, 8.1])            # V
print("R =", round(V @ I / (I @ I), 4), "ohm")  # minimos cuadrados, origen
u = np.array([1, 1, 1.0])                      # entrada constante
x = np.array([0, 2, 3.6, 4.88])                # respuesta medida
A = np.column_stack([x[:-1], u])               # renglon k: (x[k], u[k])
(a, b), *_ = np.linalg.lstsq(A, x[1:], rcond=None)
dt = 1.0                                       # s
print(f"a = {a:.4f}, b = {b:.4f}, K = {b / (1 - a):.4f},",
      f"tau = {-dt / np.log(a):.4f} s")
xs = [0.0]                                     # simulacion libre del modelo
for uk in u:
    xs.append(a * xs[-1] + b * uk)
print("simulacion:", np.round(xs, 4))
