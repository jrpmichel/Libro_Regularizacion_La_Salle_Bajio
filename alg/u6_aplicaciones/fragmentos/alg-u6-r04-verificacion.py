# ID: ALG-U6-R04
# Libro: Alg. U6 · problema resuelto ALG-U6-04 (verificacion)
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-r04-verificacion.py
import sympy as sp

R = sp.Rational
AB, BC, AD, CD, BD, Ax, Ay, Cy = sp.symbols("AB BC AD CD BD Ax Ay Cy")
c, s = R(3, 5), R(4, 5)                  # cos y sen de AD; CD es simetrica
ecs = [AB + c * AD + Ax, s * AD + Ay,                     # nodo A
       -AB + BC, BD - 12,                                 # nodo B
       -BC - c * CD, s * CD + Cy,                         # nodo C
       -c * AD + c * CD + 3, -s * AD - s * CD - BD]       # nodo D
sol = sp.solve(ecs, [AB, BC, AD, CD, BD, Ax, Ay, Cy])
print(sol)
M_A = 6 * sol[Cy] - 3 * 12 - 4 * 3      # momento respecto a A
print("Fx:", sol[Ax] + 3, " Fy:", sol[Ay] + sol[Cy] - 12, " M_A:", M_A)
A, b = sp.linear_eq_to_matrix(ecs, [AB, BC, AD, CD, BD, Ax, Ay, Cy])
print("det A =", A.det())
