# ID: ALG-U1-S4
# Libro: Alg. U1 · subtema 1.4: producto cruz, orden de los factores y area
# Repositorio: alg/u1_vectores/fragmentos/alg-u1-s4-cruz.py
import numpy as np

a = np.array([1.0, 2.0, 3.0])
b = np.array([4.0, 5.0, 6.0])
c = np.cross(a, b)
print("a x b =", c)                    # [-3.  6. -3.]
print("b x a =", np.cross(b, a))       # mismo vector, signo contrario
print("(a x b) . a =", np.dot(c, a), " (a x b) . b =", np.dot(c, b))
print("area del paralelogramo =", round(np.linalg.norm(c), 4))   # 7.3485

i, j, k = np.eye(3)
print("i x j =", np.cross(i, j), " j x i =", np.cross(j, i))   # k y -k
print("a x a =", np.cross(a, a))       # vector cero: paralelos
