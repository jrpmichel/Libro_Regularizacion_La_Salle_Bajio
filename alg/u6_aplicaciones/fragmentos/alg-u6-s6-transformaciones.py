# ID: ALG-U6-S6
# Libro: Alg. U6 · subtema 6.6: matrices homogeneas, inversa y orden
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-s6-transformaciones.py
import sympy as sp

def H(R, t):                                   # matriz homogenea 4 x 4
    T = sp.eye(4); T[:3, :3] = R; T[:3, 3] = sp.Matrix(t); return T

Rz = sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])   # 90 grados en z
T = H(Rz, [2, 0, 0])
print("punto (1,0,0):", list(T * sp.Matrix([1, 0, 0, 1]))[:3])
print("direccion (1,0,0):", list(T * sp.Matrix([1, 0, 0, 0]))[:3])
T2 = H(Rz, [1, 2, 3])
Tinv = H(Rz.T, -Rz.T * sp.Matrix([1, 2, 3]))         # [R^T, -R^T t]
print("T T^-1 = I:", T2 * Tinv == sp.eye(4))
print("T^T sirve de inversa:", T2.T * T2 == sp.eye(4))
Rx = sp.Matrix([[1, 0, 0], [0, 0, -1], [0, 1, 0]])   # 90 grados en x
Tr, Tq, p = H(sp.eye(3), [0, 0, 1]), H(Rx, [0, 0, 0]), sp.Matrix([0, 1, 0, 1])
print("traslada*gira:", list(Tr * Tq * p)[:3])
print("gira*traslada:", list(Tq * Tr * p)[:3])
