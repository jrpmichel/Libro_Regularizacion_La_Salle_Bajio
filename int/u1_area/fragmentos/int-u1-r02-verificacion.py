# ID: INT-U1-R02
# Libro: Int. U1 · problema resuelto INT-U1-02 (verificación)
# Repositorio: int/u1_area/fragmentos/int-u1-r02-verificacion.py
import numpy as np

t = np.array([0, 10, 20, 30, 40.0])      # s
P = np.array([0, 1.2, 2.0, 2.4, 2.5])    # kW
assert np.all(np.diff(P) > 0)            # lecturas crecientes
dt = np.diff(t)                          # 4 intervalos de 10 s
s = np.sum(P[:-1] * dt)                  # altura izquierda (mínimo)
S = np.sum(P[1:] * dt)                   # altura derecha (máximo)
print("s =", round(s, 6), "kJ   S =", round(S, 6), "kJ")    # 56 y 81
print("promedio:", round((s + S) / 2, 6), "kJ +/-", round((S - s) / 2, 6))
print("en kWh:", round((s + S) / 2 / 3600, 4))
