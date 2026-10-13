# ID: ALG-U4-R01
# Libro: Alg. U4 · problema resuelto ALG-U4-01 (verificacion)
# Repositorio: alg/u4_espacios/fragmentos/alg-u4-r01-verificacion.py
import sympy as sp

E = sp.Matrix([[1, 0, sp.Rational(-3, 5)], [0, 1, sp.Rational(-4, 5)]])
print("base de S0 = Nul(E):", [list(5 * w) for w in E.nullspace()])  # (3, 4, 5)
F1, F2, F3 = sp.symbols("F1:4")
print("S1:", sp.linsolve((E, sp.Matrix([0, 10])), F1, F2, F3))       # E F = -P
f, g = sp.Matrix([3, 4, 5]), sp.Matrix([-6, -8, -10])
print("suma de dos estados de S0 en S0:", E * (f + g) == sp.zeros(2, 1))
print("E (0, 0, 0) =", list(E * sp.zeros(3, 1)), "!= (0, 10): 0 no esta en S1")
F = sp.Matrix([6, 8, 10])                                           # en S2
print("2 F =", list(2 * F), " max |F_i| =", max(abs(t) for t in 2 * F), "> 10")
