# ID: ALG-U6-R06
# Libro: Alg. U6 · problema resuelto ALG-U6-06 (verificacion)
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-r06-verificacion.py
import numpy as np

u = np.array([1] * 8 + [0] * 6 + [1] * 7 + [0.5] * 4)       # u[0] ... u[24]
x = np.array([0.0, 15.8, 29.2, 40.3, 51.4, 61.1, 70.6, 78.2, 85.9, 75.8,
              70.1, 62.2, 56.5, 50.4, 45.2, 56.1, 65.8, 73.6, 81.3, 88.6,
              93.9, 99.1, 97.7, 94.8, 92.1, 90.9])          # x[0] ... x[25]
A = np.column_stack([x[:-1], u])                 # renglon k: (x[k], u[k])
(a, b), *_ = np.linalg.lstsq(A, x[1:], rcond=None)
r = x[1:] - A @ np.array([a, b])
dt = 10.0                                        # s
print(f"a = {a:.5f}, b = {b:.4f}")
print(f"tau = {-dt / np.log(a):.1f} s, K = {b / (1 - a):.1f} C")
print(f"RMSE a un paso = {np.sqrt(r @ r / r.size):.2f} C")
xs = [0.0]
for uk in u:                                     # simulacion libre
    xs.append(a * xs[-1] + b * uk)
e = x - np.array(xs)
print(f"simulacion: RMSE = {np.sqrt(e @ e / e.size):.2f} C,",
      f"maximo = {abs(e).max():.2f} C")
(av, bv), *_ = np.linalg.lstsq(A[:14], x[1:15], rcond=None)   # k = 0..13
xv = [x[14]]
for uk in u[14:]:                                # valida con k = 14..24
    xv.append(av * xv[-1] + bv * uk)
ev = x[15:] - np.array(xv[1:])
rv = np.sqrt(ev @ ev / ev.size)
print(f"validacion: a = {av:.4f}, b = {bv:.2f}, RMSE = {rv:.2f} C,",
      f"maximo = {abs(ev).max():.2f} C")
K, tau = b / (1 - a), -dt / np.log(a)
print(f"u para 120 C: {120 / K:.3f}; a plena potencia llega en",
      f"{tau * np.log(K / (K - 120)):.0f} s")
