# ID: ALG-U4-S2
# Libro: Alg. U4 · subtema 4.2: independencia con el rango y relacion de dependencia
# Repositorio: alg/u4_espacios/fragmentos/alg-u4-s2-independencia.py
import sympy as sp

v1, v2, v3 = sp.Matrix([1, 2, 0]), sp.Matrix([0, 1, 3]), sp.Matrix([2, 5, 3])
for a, b, nombre in ((v1, v2, "v1, v2"), (v1, v3, "v1, v3"), (v2, v3, "v2, v3")):
    print(nombre, "paralelos:", a.cross(b) == sp.zeros(3, 1))   # False
A = sp.Matrix.hstack(v1, v2, v3)
print("rango:", A.rank(), "de", A.cols)          # 2 de 3: dependientes
c = A.nullspace()[0]
c = c / sp.gcd(list(c))                          # coeficientes enteros
c = -c if c[0] < 0 else c                        # primer coeficiente positivo
print("relacion c1 v1 + c2 v2 + c3 v3 = 0 con c =", list(c))
print("comprobacion:", list(A * c))              # [0, 0, 0]
print("det:", A.det())                           # 0
