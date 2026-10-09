# ID: INT-U1-R03
# Libro: Int. U1 · problema resuelto INT-U1-03 (verificación)
# Repositorio: int/u1_area/fragmentos/int-u1-r03-verificacion.py
import numpy as np

tv = np.array([0, 2, 8, 10, 12, 13, 16, 17.0])     # s
vv = np.array([0, 1.2, 1.2, 0, 0, -0.8, -0.8, 0])  # m/s
t = np.linspace(0, 17, 17001)                      # incluye los vértices
v = np.interp(t, tv, vv)
pasos = (v[1:] + v[:-1]) / 2 * np.diff(t)          # trapecios delgados
s = np.concatenate([[0], np.cumsum(pasos)])        # posición acumulada
print("posición final:", round(s[-1], 6), "m")                # 6.4
print("odómetro:", round(np.sum(np.abs(pasos)), 6), "m")      # 12.8
print("posición máxima:", round(s.max(), 6), "m")             # 9.6
assert np.allclose([s[-1], np.abs(pasos).sum(), s.max()], [6.4, 12.8, 9.6])
