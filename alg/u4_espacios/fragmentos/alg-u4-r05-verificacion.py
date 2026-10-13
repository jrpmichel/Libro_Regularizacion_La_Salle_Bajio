# ID: ALG-U4-R05
# Libro: Alg. U4 · problema resuelto ALG-U4-05 (verificacion)
# Repositorio: alg/u4_espacios/fragmentos/alg-u4-r05-verificacion.py
import sympy as sp

R = sp.Rational
P1 = sp.Matrix([[1, 0, 0], [0, 1, 0]])                 # camara 1
c, s = R(4, 5), R(3, 5)
Ry = sp.Matrix([[c, 0, s], [0, 1, 0], [-s, 0, c]])     # R_y(phi), mano derecha
P2 = P1 * Ry.T                                         # camara 2: P1 R_y(-phi)
print("P2 =", P2.tolist())                             # [[4/5, 0, -3/5], [0, 1, 0]]
print("nucleos:", list(P1.nullspace()[0]), list(4 * P2.nullspace()[0]))   # vision
Q, Qp = sp.Matrix([2, 1, 5]), sp.Matrix([2, 1, -1])
print("camara 1:", list(P1 * Q), list(P1 * Qp))        # iguales
print("camara 2:", list(P2 * Q), list(P2 * Qp))        # distintas
M = P1.col_join(P2)                                    # las dos camaras juntas
print("rango de M:", M.rank(), " nucleo:", M.nullspace())   # 3 y []
m = sp.Matrix([2, 1, R(-7, 5), 1])                     # imagenes medidas
print("reconstruccion:", sp.linsolve((M, m), sp.symbols("x y z")))  # (2, 1, 5)
Pz = sp.Matrix([[c, s, 0], [-s, c, 0]])                # girada sobre z
print("rango con giro sobre z:", P1.col_join(Pz).rank())   # 2
