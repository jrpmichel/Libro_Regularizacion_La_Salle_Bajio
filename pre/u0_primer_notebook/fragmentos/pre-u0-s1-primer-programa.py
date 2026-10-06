# ID: PRE-U0-S1
# Libro: Preliminares · capítulo Python y Google Colab: primer programa
# Repositorio: pre/u0_primer_notebook/fragmentos/pre-u0-s1-primer-programa.py
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x**2 + 1          # la regla; cámbiala y vuelve a ejecutar

x = np.arange(0, 5)          # entradas: 0, 1, 2, 3, 4
y = f(x)                     # salidas

print("x =", x)
print("y =", y)

plt.plot(x, y, "o-")
plt.grid(True)
plt.xlabel("x")
plt.ylabel("y")
plt.show()
