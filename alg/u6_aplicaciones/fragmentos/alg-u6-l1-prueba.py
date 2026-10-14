# ID: ALG-U6-L1
# Libro: Alg. U6 · Laboratorio del Error: el momento al reves
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-l1-prueba.py
import numpy as np

r = np.array([0.25, 0.0, 0.10])     # m, del centro de la muneca a la punta
F = np.array([0.0, 40.0, -30.0])    # N, fuerza de la pieza sobre la punta
bien, mal = np.cross(r, F), np.cross(F, r)
print("r x F =", bien, " F x r =", mal)
seno = np.sqrt(1 - (r @ F / (np.linalg.norm(r) * np.linalg.norm(F)))**2)
print("|r||F| sen(theta) =",
      round(np.linalg.norm(r) * np.linalg.norm(F) * seno, 4))
print("normas:", round(np.linalg.norm(bien), 4),
      round(np.linalg.norm(mal), 4))
print("termino x por y:", r[0] * F[1], "N m en +z (i x j = k)")
for nombre, M in (("bien", bien), ("al reves", mal)):
    motor = -M[2]                        # par que compensa el momento en z
    print(nombre, ": par del motor", motor, "N m; par neto",
          bien[2] + motor, "N m")
