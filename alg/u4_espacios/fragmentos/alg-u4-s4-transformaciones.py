# ID: ALG-U4-S4
# Libro: Alg. U4 · subtema 4.4: rotacion, proyeccion y reflexion como matrices
# Repositorio: alg/u4_espacios/fragmentos/alg-u4-s4-transformaciones.py
import numpy as np
import sympy as sp

def R(grados):                           # rotacion antihoraria en el plano
    t = np.deg2rad(grados)               # grados -> radianes
    return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])

print("R(30) (2, 0) =", np.round(R(30) @ [2, 0], 4))      # (1.7321, 1.0)
print("R(25) R(35) = R(60):", np.allclose(R(25) @ R(35), R(60)))
u = sp.Matrix([3, 4])
P = u * u.T / u.dot(u)                   # proyeccion sobre la recta de u
H = 2 * P - sp.eye(2)                    # reflexion respecto a esa recta
print("P =", P.tolist(), " P^2 = P:", P * P == P)
print("P (5, 0) =", list(P * sp.Matrix([5, 0])), " P (4, -3) =",
      list(P * sp.Matrix([4, -3])))
w = P.nullspace()[0]
print("nucleo: gen", list(-3 * w), " rango:", P.rank())     # gen{(4, -3)}
print("H (5, 0) =", list(H * sp.Matrix([5, 0])), " det H =", H.det())
