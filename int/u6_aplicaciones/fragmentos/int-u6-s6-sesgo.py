# ID: INT-U6-S6
# Libro: Int. U6 · subtema 6.6: doble integración con y sin sesgo
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-s6-sesgo.py
import numpy as np

dt = 0.01
t = np.arange(0, 60 + dt / 2, dt)                   # s
tm = t[:-1] + dt / 2                                # puntos medios

def a_real(t):
    """Acelera 0.5 m/s^2 hasta 2 s, frena de 18 a 20 s."""
    return np.where(t < 2, 0.5, 0.0) - np.where((t > 18) & (t < 20), 0.5, 0.0)

def posicion(sesgo):
    v = np.concatenate(([0.0], np.cumsum((a_real(tm) + sesgo) * dt)))
    return np.concatenate(([0.0], np.cumsum((v[1:] + v[:-1]) / 2 * dt)))

s_real, s_medida = posicion(0.0), posicion(0.05)    # sesgo de 0.05 m/s^2
for k in (2000, 6000):                              # t = 20 s y 60 s
    print(f"t = {t[k]:.0f} s: real {s_real[k]:.2f} m,"
          f" medida {s_medida[k]:.2f} m")
