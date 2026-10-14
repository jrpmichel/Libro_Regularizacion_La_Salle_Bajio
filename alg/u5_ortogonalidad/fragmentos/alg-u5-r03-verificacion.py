# ID: ALG-U5-R03
# Libro: Alg. U5 · problema resuelto ALG-U5-03 (verificacion)
# Repositorio: alg/u5_ortogonalidad/fragmentos/alg-u5-r03-verificacion.py
import numpy as np

Rt = np.array([[0.80, -0.59, 0.01],      # columnas: ejes x, y, z del dron
               [0.60, 0.81, -0.02],
               [0.00, 0.02, 1.00]])
print("Rt^T Rt - I =\n", np.round(Rt.T @ Rt - np.eye(3), 5))

def gs(M, orden):                         # Gram-Schmidt en el orden dado
    Q, hechos = np.zeros_like(M), []
    for j in orden:
        u = M[:, j] - sum((M[:, j] @ q) * q for q in hechos)
        hechos.append(u / np.linalg.norm(u))
        Q[:, j] = hechos[-1]
    return Q

for orden, nombre in (([0, 1, 2], "x, y, z"), ([2, 0, 1], "z, x, y")):
    Q = gs(Rt, orden)
    giro = [np.degrees(np.arccos(Q[:, j] @ Rt[:, j] / np.linalg.norm(Rt[:, j])))
            for j in range(3)]
    error = abs(Q.T @ Q - np.eye(3)).max()
    print(f"orden {nombre}: giros x, y, z = {np.round(giro, 2)} grados")
    print(f"   det = {np.linalg.det(Q):.6f}, max|Q^T Q - I| = {error:.0e}")
