# ID: DIF-U1-S7
# Libro: Dif. U1 · subtema 1.6: seno y coseno
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-s7-seno-y-coseno.py
import numpy as np

th = np.pi / 6                                    # 30 grados, expresados en radianes
print(round(np.cos(th), 4), round(np.sin(th), 4)) # 0.866 0.5: punto del círculo unitario
print(round(np.sin(30), 4), round(np.sin(np.deg2rad(30)), 4))   # -0.988 (30 rad) frente a 0.5 (30 grados)
t = np.linspace(0, 2*np.pi, 9)
print(np.round(np.sin(t)**2 + np.cos(t)**2, 12))  # identidad pitagórica: todos 1
A, w, fi = 8, 100*np.pi, np.pi/4                  # v(t) = A sen(w t + fi)
print("periodo =", 2*np.pi/w, "s | frecuencia =", w/(2*np.pi), "Hz")   # 0.02 s | 50 Hz
