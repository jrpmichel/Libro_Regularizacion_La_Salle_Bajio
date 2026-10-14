# ID: ALG-U6-S3
# Libro: Alg. U6 · subtema 6.3: mezclas posibles con tres chatarras
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-s3-mezclas.py
import numpy as np
import matplotlib.pyplot as plt

S = np.array([[1.00, 1.00, 1.00],        # total
              [0.80, 0.60, 0.65],        # fraccion de cobre de cada chatarra
              [0.18, 0.35, 0.25]])       # fraccion de zinc
cu, zn = np.meshgrid(np.linspace(55, 85, 61), np.linspace(14, 40, 53))
si = []
for c, z in zip(cu.ravel(), zn.ravel()):
    x = np.linalg.solve(S, [1000, 10 * c, 10 * z])   # kg para 1000 kg
    si.append(np.all(x >= -1e-9))                    # factible: x >= 0
si = np.array(si)
print("objetivos factibles:", si.sum(), "de", si.size)
plt.plot(cu.ravel()[si], zn.ravel()[si], "k.", ms=3, label="factible")
plt.plot(cu.ravel()[~si], zn.ravel()[~si], "x", c="0.7", ms=3,
         label="alguna cantidad negativa")
plt.xlabel("cobre (%)"); plt.ylabel("zinc (%)"); plt.legend(); plt.show()
