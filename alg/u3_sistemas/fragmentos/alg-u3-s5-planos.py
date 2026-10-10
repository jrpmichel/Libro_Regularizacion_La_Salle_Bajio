# ID: ALG-U3-S5
# Libro: Alg. U3 · subtema 3.5: tres planos y su punto comun
# Repositorio: alg/u3_sistemas/fragmentos/alg-u3-s5-planos.py
import numpy as np
import matplotlib.pyplot as plt

A = np.array([[1.0, 1, 1], [1, -1, 2], [2, 1, -1]])
b = np.array([3.0, 2, 2])
p = np.linalg.solve(A, b)                       # punto comun
print("punto comun:", p)                        # [1. 1. 1.]
u, v = np.meshgrid(np.linspace(-1, 3, 9), np.linspace(-1, 3, 9))
ax = plt.figure().add_subplot(projection="3d")
for fila, bi, estilo in zip(A, b, ["-", "--", ":"]):
    w = (bi - fila[0] * u - fila[1] * v) / fila[2]   # despeja z
    ax.plot_wireframe(u, v, w, linestyle=estilo, linewidth=0.6,
                      color="k")
ax.scatter(*p, s=60, color="k")
ax.set_xlabel("x")
ax.set_ylabel("y")
ax.set_zlabel("z")
plt.show()
