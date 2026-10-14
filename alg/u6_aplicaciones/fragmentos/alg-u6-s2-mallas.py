# ID: ALG-U6-S2
# Libro: Alg. U6 · subtema 6.2: dos mallas, signos de las fuentes y comprobacion
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-s2-mallas.py
import sympy as sp

R1, R3, R2 = 2, 4, 4                     # ohm; R3 es la rama compartida
Z = sp.Matrix([[R1 + R3, -R3],
               [-R3, R2 + R3]])          # matriz de resistencias de malla
for V2 in (0, 4):                    # V2 con + arriba: resta en la malla 2
    V = sp.Matrix([10, -V2])
    i = Z.solve(V)
    print(f"V2 = {V2} V: i1 = {i[0]} A, i2 = {i[1]} A,",
          f"por R3 bajan {i[0] - i[1]} A")
Req = R1 + sp.Rational(R3 * R2, R3 + R2)  # serie-paralelo con V2 = 0
print("comprobacion: 10 V /", Req, "ohm =", sp.Rational(10) / Req, "A")
