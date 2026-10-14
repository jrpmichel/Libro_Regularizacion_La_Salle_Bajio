# ID: PRE-U0-S2
# Libro: Anexo C · numpy: arreglos, operaciones, algebra lineal y angulos
# Repositorio: pre/u0_primer_notebook/fragmentos/pre-u0-s2-numpy.py
import numpy as np

x = np.linspace(0, 2, 5)              # 5 puntos de 0 a 2, extremos incluidos
y = x**2 + 1                          # se aplica a cada elemento
print("x =", x, " y =", y)
a, b = np.array([1, 2, 2]), np.array([2, 1, -2])
print("a . b =", a @ b, "  a x b =", np.cross(a, b))
print("|a| =", np.linalg.norm(a))
A = np.array([[2.0, 1.0], [1.0, 3.0]])
print("solucion de A x = (3, 5):", np.linalg.solve(A, [3, 5]))
print("det A =", round(np.linalg.det(A), 10))
print("sen 30 grados =", np.sin(np.radians(30)), "  np.sin(30) =", np.sin(30))
print("0.1 + 0.2 == 0.3:", 0.1 + 0.2 == 0.3,
      "  np.isclose:", np.isclose(0.1 + 0.2, 0.3))
