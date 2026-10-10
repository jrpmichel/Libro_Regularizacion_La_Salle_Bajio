# ID: ALG-U2-S4
# Libro: Alg. U2 · subtema 2.4: imagen del cuadrado unitario y area = |det A|
# Repositorio: alg/u2_matrices/fragmentos/alg-u2-s4-geometria.py
import numpy as np
import matplotlib.pyplot as plt

A = np.array([[2.0, 1.0], [1.0, 1.0]])
cuadrado = np.array([[0, 1, 1, 0], [0, 0, 1, 1]])   # vertices como columnas
imagen = A @ cuadrado

x, y = imagen
# area de un poligono a partir de sus vertices (formula del zapato)
area = 0.5 * abs(np.dot(x, np.roll(y, -1)) - np.dot(y, np.roll(x, -1)))
print("vertices de la imagen:\n", imagen)
print("area de la imagen =", area, " det A =", round(np.linalg.det(A), 10))

plt.fill(*cuadrado, fill=False, ls="--", label="cuadrado unitario")
plt.fill(*imagen, alpha=0.3, label="imagen bajo A")
plt.axis("equal"); plt.grid(True); plt.legend(); plt.show()
