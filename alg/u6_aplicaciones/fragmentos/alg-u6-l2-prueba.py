# ID: ALG-U6-L2
# Libro: Alg. U6 · Laboratorio del Error: el balance que siempre cierra
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-l2-prueba.py
import numpy as np

R1, R2, R3 = 4.0, 3.0, 2.0               # ohm; R3 compartida
Z = np.array([[R1 + R3, -R3], [-R3, R2 + R3]])
for nombre, V in (("signo equivocado", [12, 6]), ("correcto", [12, -6])):
    i = np.linalg.solve(Z, V)
    pR = R1 * i[0]**2 + R2 * i[1]**2 + R3 * (i[0] - i[1])**2
    print(f"{nombre}: i = {np.round(i, 4)} A,",
          f"P_R = {pR:.2f} W, P_fuentes = {np.dot(V, i):.2f} W")
rng = np.random.default_rng(1)           # cualquier lado derecho cierra
for _ in range(3):
    V = rng.uniform(-20, 20, 2); i = np.linalg.solve(Z, V)
    print("V =", np.round(V, 2), " i^T Z i - i^T V =",
          round(i @ Z @ i - i @ V, 12))
i_sola = np.linalg.solve(Z, [0, -6])     # solo V2: debe empujar al reves
print("solo V2:", np.round(i_sola, 4), "A (i2 < 0: antihorario)")
