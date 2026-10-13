# ID: ALG-U4-S3
# Libro: Alg. U4 · subtema 4.3: coordenadas, cambio de base y teorema de la dimension
# Repositorio: alg/u4_espacios/fragmentos/alg-u4-s3-coordenadas.py
import sympy as sp

PB = sp.Matrix([[1, 1], [0, 1]])        # columnas: b1 = (1, 0), b2 = (1, 1)
x = sp.Matrix([3, 2])
cB = PB.solve(x)                         # [x]_B = PB^(-1) x
print("[x]_B =", list(cB), " regreso:", list(PB * cB))
PC = sp.Matrix([[1, 1], [1, -1]])       # otra base C
PCB = PC.inv() * PB                      # de coordenadas B a coordenadas C
print("P_C<-B =", PCB.tolist(), " [x]_C =", list(PCB * cB))

A = sp.Matrix([[1, 2, 3], [2, 4, 7], [3, 6, 10]])
R, piv = A.rref()
print("base de Col(A):", [list(A[:, j]) for j in piv])    # columnas de A
print("base de Nul(A):", [list(w) for w in A.nullspace()])
print("rango + dim Nul =", len(piv), "+", len(A.nullspace()), "=", A.cols)
