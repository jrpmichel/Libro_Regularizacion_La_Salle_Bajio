# ID: INT-U2-R02
# Libro: Int. U2 · problema resuelto INT-U2-02 (verificación)
# Repositorio: int/u2_riemann/fragmentos/int-u2-r02-verificacion.py
import sympy as sp

i, n = sp.symbols("i n", integer=True, positive=True)
k, d = 2000, sp.Rational(15, 100)              # N/m y m
dx = d / n
R_n = sp.simplify(sp.summation(k * (i * dx) * dx, (i, 1, n)))
print("R_n =", sp.factor(R_n))                  # 45(n + 1)/(2n)
print("W =", sp.limit(R_n, n, sp.oo), "J")      # 45/2
L10 = sp.summation(k * (i * d / 10) * (d / 10), (i, 0, 9))
print("R_10 =", float(R_n.subs(n, 10)), "J   L_10 =", float(L10), "J")
