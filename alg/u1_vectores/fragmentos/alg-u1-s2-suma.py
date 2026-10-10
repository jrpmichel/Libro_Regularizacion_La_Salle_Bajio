# ID: ALG-U1-S2
# Libro: Alg. U1 · subtema 1.2: suma de dos fuerzas y vector unitario de la resultante
# Repositorio: alg/u1_vectores/fragmentos/alg-u1-s2-suma.py
import numpy as np
import matplotlib.pyplot as plt

F1 = np.array([300.0, 0.0])       # N, cuerda 1
F2 = np.array([0.0, 400.0])       # N, cuerda 2, perpendicular a la 1
R = F1 + F2                       # se suman componente a componente
print("R =", R, "N")
print("|R| =", np.linalg.norm(R), "N")                       # 500, no 700
print("|F1| + |F2| =", np.linalg.norm(F1) + np.linalg.norm(F2), "N")
u = R / np.linalg.norm(R)         # unitario: misma direccion, magnitud 1
print("u =", u, " |u| =", np.linalg.norm(u))                 # [0.6 0.8], 1.0

o = np.zeros(2)
for v, c, e in ((F1, "gray", "F1"), (F2, "gray", "F2"), (R, "navy", "R")):
    plt.quiver(*o, *v, angles="xy", scale_units="xy", scale=1,
               color=c, label=e)
plt.plot([F1[0], R[0]], [F1[1], R[1]], "k:")   # paralelogramo
plt.plot([F2[0], R[0]], [F2[1], R[1]], "k:")
plt.axis("equal"); plt.xlim(-50, 450); plt.ylim(-50, 450)
plt.legend(); plt.grid(True)
plt.xlabel("x (N)"); plt.ylabel("y (N)"); plt.show()
