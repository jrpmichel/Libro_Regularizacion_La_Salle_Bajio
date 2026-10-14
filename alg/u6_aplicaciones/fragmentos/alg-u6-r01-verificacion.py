# ID: ALG-U6-R01
# Libro: Alg. U6 · problema resuelto ALG-U6-01 (verificacion)
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-r01-verificacion.py
import sympy as sp

R = sp.Rational
P = sp.Matrix([0, 0, 0])                                  # anillo (m)
anc = {"A": sp.Matrix([R(-18, 10), 0, R(24, 10)]),
       "B": sp.Matrix([R(108, 100), R(144, 100), R(24, 10)]),
       "D": sp.Matrix([0, -1, R(24, 10)])}
u = {k: (Q - P) / (Q - P).norm() for k, Q in anc.items()}
print("longitudes (m):", [(Q - P).norm() for Q in anc.values()])
U = sp.Matrix.hstack(u["A"], u["B"], u["D"])
T = U.solve(sp.Matrix([0, 0, R(38, 10)]))                 # U T = (0, 0, W)
print("T_A, T_B, T_D (kN):", list(T), [float(t) for t in T])
F = -T[2] * u["D"]                       # la cadena jala la mensula hacia P
r = anc["D"] - sp.Matrix([0, R(-16, 10), R(24, 10)])      # de C a D
print("F en D (kN):", list(F), " r (m):", list(r))
print("M_C = r x F (kN m):", list(r.cross(F)))
print("suma de fuerzas en P:", list(U * T - sp.Matrix([0, 0, R(38, 10)])))
