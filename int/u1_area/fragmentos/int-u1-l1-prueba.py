# ID: INT-U1-L1
# Libro: Int. U1 · Laboratorio del Error 1.1 (prueba con código)
# Repositorio: int/u1_area/fragmentos/int-u1-l1-prueba.py
import numpy as np

tv = np.array([0, 1, 9, 10, 15, 16, 21, 22.0])   # s; pausa de 10 a 15 s
vv = np.array([0, 2, 2, 0, 0, -2, -2, 0])        # m/s
t = np.linspace(0, 22, 22001)
v = np.interp(t, tv, vv)
pasos = (v[1:] + v[:-1]) / 2 * np.diff(t)
print("área con signo:", round(pasos.sum(), 6), "m")          # 6: desplazamiento
print("área total:", round(np.abs(pasos).sum(), 6), "m")      # 30: recorrido del cable
