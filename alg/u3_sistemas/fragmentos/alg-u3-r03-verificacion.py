# ID: ALG-U3-R03
# Libro: Alg. U3 · problema resuelto ALG-U3-03 (verificacion)
# Repositorio: alg/u3_sistemas/fragmentos/alg-u3-r03-verificacion.py
import numpy as np

A = np.array([[4, -1, -1, 0], [-1, 4, 0, -1],
              [-1, 0, 4, -1], [0, -1, -1, 4]], float)
b = np.array([140, 140, 60, 60], float)
T = np.array([60, 60, 40, 40], float)               # reportado
print("residuo b - A T:", b - A @ T)                # [0. 0. 0. 0.]
print("det A =", round(np.linalg.det(A), 6))        # 192
print("det forma promedio =", round(np.linalg.det(A / 4), 6))   # 0.75

def gauss(n):
    return (n**3 + 3 * n**2 - n) // 3

def cofactores(n):
    return 2 if n == 2 else n * (cofactores(n - 1) + 1)

print("n = 4: Cramer", 5 * cofactores(4) + 4, "Gauss", gauss(4))
print("n = 100: Gauss", gauss(100))                 # 343300
N = 10                                              # malla 10 x 10
L = 4 * np.eye(N * N)
for k in range(N * N):
    i, j = divmod(k, N)
    for ii, jj in ((i + 1, j), (i - 1, j), (i, j + 1), (i, j - 1)):
        if 0 <= ii < N and 0 <= jj < N:
            L[k, ii * N + jj] = -1
print("malla 10x10, det promedio:", np.linalg.det(L / 4))   # 2.5e-8
print("numero de condicion:", round(np.linalg.cond(L), 1))  # 48.4
