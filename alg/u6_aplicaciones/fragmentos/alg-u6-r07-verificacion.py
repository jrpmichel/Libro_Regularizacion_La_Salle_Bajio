# ID: ALG-U6-R07
# Libro: Alg. U6 · problema resuelto ALG-U6-07 (verificacion)
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-r07-verificacion.py
import numpy as np

R = np.column_stack([[0, -1, 0], [0, 0, -1], [1, 0, 0]])   # x_C, y_C, z_C
t = np.array([-0.3, 0.0, 0.5])                             # m
print("R^T R = I:", np.allclose(R.T @ R, np.eye(3)),
      " det =", round(np.linalg.det(R), 12))
T = np.eye(4); T[:3, :3] = R; T[:3, 3] = t
pC = np.array([0.1, 0.2, 0.9, 1.0])
print("p_B =", (T @ pC)[:3])
print("con R^T:", R.T @ pC[:3] + t, "  R(p + t):", R @ (pC[:3] + t))
Ti = np.eye(4); Ti[:3, :3] = R.T; Ti[:3, 3] = -R.T @ t
print("base vista desde la camara:", Ti[:3, 3], " T Ti = I:",
      np.allclose(T @ Ti, np.eye(4)))
RF = np.diag([1, -1, -1])           # brida con la herramienta hacia abajo
pF = (T @ pC)[:3] - RF @ np.array([0, 0, 0.15])
print("brida:", pF)
