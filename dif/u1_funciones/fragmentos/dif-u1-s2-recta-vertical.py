# ID: DIF-U1-S2
# Libro: Dif. U1 · subtema 1.2: criterio de la recta vertical
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-s2-recta-vertical.py
import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0, 2*np.pi, 400)
x = np.linspace(-5, 5, 400)
plt.plot(5*np.cos(t), 5*np.sin(t), label="x^2 + y^2 = 25")
plt.plot(x, np.sqrt(25 - x**2), "--", label="y = sqrt(25 - x^2)")
plt.axvline(3, color="gray")              # recta vertical x = 3
plt.axis("equal"); plt.legend(); plt.grid(True); plt.show()
c = 3                                     # cuántas y le tocan a x = c
print("y en x = 3:", np.sqrt(25 - c**2), -np.sqrt(25 - c**2))   # 4.0 -4.0
