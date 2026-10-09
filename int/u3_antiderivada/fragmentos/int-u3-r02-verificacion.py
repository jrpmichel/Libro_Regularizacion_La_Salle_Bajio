# ID: INT-U3-R02
# Libro: Int. U3 · problema resuelto INT-U3-02 (verificación)
# Repositorio: int/u3_antiderivada/fragmentos/int-u3-r02-verificacion.py
import sympy as sp

q, K = sp.symbols("q K")
Cm = 120 - sp.Rational(2, 5) * q          # costo marginal, pesos/pieza
C = sp.integrate(Cm, q) + K               # costo total, pesos
print("C(q) =", C)
CA, CB = C.subs(K, 8000), C.subs(K, 11000)
print("C_B - C_A =", sp.simplify(CB - CA))                    # 3000
print("C(150) - C(100) =", C.subs(q, 150) - C.subs(q, 100))   # 3500
print("C'(100), C'(150) =", Cm.subs(q, 100), Cm.subs(q, 150)) # 80, 60
