# ID: ALG-U2-S1
# Libro: Alg. U2 · subtema 2.1: producto de matrices, orden de los factores y transpuesta
# Repositorio: alg/u2_matrices/fragmentos/alg-u2-s1-producto.py
import numpy as np

A = np.array([[2, 0, 1],
              [1, -1, 3]])          # 2x3
B = np.array([[1, 2],
              [0, 1],
              [4, -2]])             # 3x2
print("A es", A.shape, " B es", B.shape)
print("AB =\n", A @ B)              # 2x2: renglon de A por columna de B
print("BA =\n", B @ A)              # 3x3: otra matriz, de otro tamano
print("(AB)^T == B^T A^T:", np.array_equal((A @ B).T, B.T @ A.T))

try:
    A @ A                           # 2x3 por 2x3: no esta definido
except ValueError:
    print("A A no existe: 2x3 por 2x3, las columnas de A no son")
    print("tantas como los renglones de A")
