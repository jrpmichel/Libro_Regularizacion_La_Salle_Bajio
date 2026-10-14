# ID: ALG-U6-S4
# Libro: Alg. U6 · subtema 6.4: metodo de nodos para la armadura de la semilla
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-s4-armadura.py
import numpy as np

nodos = {"A": (0, 0), "B": (3, 0), "C": (6, 0), "D": (3, 4)}     # m
barras = [("A", "B"), ("B", "C"), ("A", "D"), ("C", "D"), ("B", "D")]
apoyos = [("A", 0), ("A", 1), ("C", 1)]          # Ax, Ay, Cy (0 = x, 1 = y)
idx = {k: i for i, k in enumerate(nodos)}
A = np.zeros((8, 8)); b = np.zeros(8)
for j, (p, q) in enumerate(barras):              # tension positiva
    d = np.subtract(nodos[q], nodos[p]); u = d / np.hypot(*d)
    A[2 * idx[p]:2 * idx[p] + 2, j] += u         # la barra jala a p hacia q
    A[2 * idx[q]:2 * idx[q] + 2, j] -= u         # y a q hacia p
for k, (n, c) in enumerate(apoyos):
    A[2 * idx[n] + c, 5 + k] = 1
for escala in (1, 2):                            # 12 kN en B, luego 24 kN
    b[2 * idx["B"] + 1] = 12 * escala            # sum F = 0: A x = -carga
    x = np.linalg.solve(A, b)
    nombres = ["AB", "BC", "AD", "CD", "BD", "Ax", "Ay", "Cy"]
    print(escala * 12, "kN:",
          ", ".join(f"{n} = {v:g}" for n, v in zip(nombres, x)))
