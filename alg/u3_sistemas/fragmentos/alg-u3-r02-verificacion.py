# ID: ALG-U3-R02
# Libro: Alg. U3 · problema resuelto ALG-U3-02 (verificacion)
# Repositorio: alg/u3_sistemas/fragmentos/alg-u3-r02-verificacion.py
import numpy as np

W = 500.0                                           # N
al, be = np.radians(20), np.radians(30)
A = np.array([[-np.cos(al), np.cos(be)], [np.sin(al), np.sin(be)]])
b = np.array([0.0, W])
d = np.linalg.det(A)
A1, A2 = A.copy(), A.copy()
A1[:, 0], A2[:, 1] = b, b                           # Cramer
T1, T2 = np.linalg.det(A1) / d, np.linalg.det(A2) / d
print(f"det A = {d:.4f}")                           # -0.7660
print(f"T1 = {T1:.1f} N, T2 = {T2:.1f} N")          # 565.3, 613.3
print("A T = b:", np.allclose(A @ [T1, T2], b))      # True
t = np.radians(5)
print(f"cables a 5 grados: {W * np.cos(t) / np.sin(2 * t):.0f} N")
