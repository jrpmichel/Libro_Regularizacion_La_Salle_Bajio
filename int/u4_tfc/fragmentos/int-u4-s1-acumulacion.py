# ID: INT-U4-S1
# Libro: Int. U4 · subtema 4.1: función de acumulación y su derivada
# Repositorio: int/u4_tfc/fragmentos/int-u4-s1-acumulacion.py
import numpy as np
import matplotlib.pyplot as plt

f = lambda t: 2 + np.cos(t)
t = np.linspace(0, 4, 2001)
dt = t[1] - t[0]
# G(x) = área de 0 a x, acumulada con trapecios
G = np.concatenate([[0], np.cumsum((f(t[:-1]) + f(t[1:])) / 2 * dt)])
dG = np.gradient(G, dt)                   # derivada numérica de G
error = np.max(np.abs(dG[1:-1] - f(t[1:-1])))   # puntos interiores
print("máx |G' - f| =", error)            # del orden de 1e-6
print("G(2.5) =", round(np.interp(2.5, t, G), 4), " exacto:",
      round(5 + np.sin(2.5), 4))          # G(x) = 2x + sen x

plt.plot(t, f(t), label="f(t) = 2 + cos t")
plt.plot(t, G, label="G(x)")
plt.plot(t, dG, ":", label="G'(x) numérica")
plt.xlabel("x"); plt.legend(); plt.grid(True); plt.show()
