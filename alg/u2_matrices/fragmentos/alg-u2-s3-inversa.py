# ID: ALG-U2-S3
# Libro: Alg. U2 · subtema 2.3: adjunta, inversa y comprobacion A A^-1 = I
# Repositorio: alg/u2_matrices/fragmentos/alg-u2-s3-inversa.py
import sympy as sp

A = sp.Matrix([[3, 1], [5, 2]])
print("det A =", A.det())                       # 1
print("adj A =", A.adjugate())                  # [[2, -1], [-5, 3]]
Ainv = A.adjugate() / A.det()
print("A^-1 =", Ainv)
print("A A^-1 = I:", A * Ainv == sp.eye(2))

B = sp.Matrix([[1, 2, 0], [0, 1, 1], [2, 0, 1]])
print("det B =", B.det(), " det(2B) =", (2 * B).det())   # 5 y 40 = 2**3 * 5
C = sp.Matrix([[2, 1, 0], [0, 2, 1], [1, 0, 2]])
print("det(BC) = det B det C:", (B * C).det() == B.det() * C.det())   # 45
S = sp.Matrix([[1, 2], [2, 4]])
print("det S =", S.det(), "-> sin inversa")     # renglones proporcionales
