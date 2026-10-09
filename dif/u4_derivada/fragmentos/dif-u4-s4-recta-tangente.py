# ID: DIF-U4-S4
# Libro: Dif. U4 · subtema 4.4: recta tangente a 1/x en x = 2
# Repositorio: dif/u4_derivada/fragmentos/dif-u4-s4-recta-tangente.py
import numpy as np
import matplotlib.pyplot as plt

f = lambda x: 1 / x
a, m = 2, -1 / 4                 # f'(2) = -1/2**2 (sección 4.3)
tangente = lambda x: f(a) + m * (x - a)

for xx in [2.1, 2.5, 3]:         # cerca de a, la recta casi es la curva
    print(xx, round(f(xx), 4), round(tangente(xx), 4))

x = np.linspace(0.6, 4, 200)
plt.plot(x, f(x), label="y = 1/x")
plt.plot(x, tangente(x), "--", label="tangente en x = 2")
plt.plot(a, f(a), "o"); plt.grid(True); plt.legend(); plt.show()
