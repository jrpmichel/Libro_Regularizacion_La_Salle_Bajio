# ID: ALG-U5-S4
# Libro: Alg. U5 · subtema 5.4: recta de minimos cuadrados del resorte y sus residuos
# Repositorio: alg/u5_ortogonalidad/fragmentos/alg-u5-s4-minimos.py
import numpy as np
import matplotlib.pyplot as plt

m = np.array([0, 100, 200, 300, 400])           # masa colgada (g)
y = np.array([10.0, 12.1, 13.8, 16.2, 17.9])    # longitud del resorte (cm)
A = np.column_stack([np.ones_like(m), m])
b0, b1 = np.linalg.solve(A.T @ A, A.T @ y)      # ecuaciones normales
r = y - (b0 + b1 * m)
print(f"y = {b0:.2f} + {b1:.4f} m;  residuos {np.round(r, 2)}")
print(f"suma de residuos {r.sum():.1e}")
print("pasa por los promedios:", np.isclose(b0 + b1 * m.mean(), y.mean()))
alarga = b1                 # 100 g alargan 100*b1 cm, que son b1 metros
print("k = 0.981 N /", round(alarga, 4), "m =", round(0.981 / alarga, 1), "N/m")
plt.plot(m, y, "ko", label="datos")
plt.plot(m, b0 + b1 * m, "k--", label="minimos cuadrados")
plt.vlines(m, y, b0 + b1 * m, colors="gray")    # residuos
plt.xlabel("Masa colgada (g)"); plt.ylabel("Longitud (cm)")
plt.legend(); plt.show()
