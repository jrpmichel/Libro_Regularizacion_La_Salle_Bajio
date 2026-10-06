# ID: DIF-U1-S5
# Libro: Dif. U1 · subtema 1.5: recta y cuadrática
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-s5-recta-y-parabola.py
import numpy as np

a, b, c = 1, -4, 3                          # f(x) = a x^2 + b x + c, con a distinto de 0
disc = b**2 - 4*a*c                         # discriminante
xv = -b / (2*a); yv = a*xv**2 + b*xv + c    # vértice
print("discriminante:", disc, "| vértice:", (xv, yv))
if disc >= 0:
    print("raíces:", (-b + np.array([1, -1])*np.sqrt(disc)) / (2*a))   # ambas, con su signo
else:
    print("sin raíces reales")
