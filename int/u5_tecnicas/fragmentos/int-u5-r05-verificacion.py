# ID: INT-U5-R05
# Libro: Int. U5 · problema resuelto INT-U5-05 (verificación)
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-r05-verificacion.py
import sympy as sp

P = sp.symbols("P", positive=True)
k = sp.Rational(3, 5)                          # 0.6 por mes
print(sp.apart(1 / (P * (1 - P)), P))          # -1/(P - 1) + 1/P
T = sp.integrate(1 / (k * P * (1 - P)),
                 (P, sp.Rational(1, 10), sp.Rational(9, 10)))
print("T =", sp.simplify(T), "=", round(float(T), 2), "meses")
# T = 2 ln 9 / 0.6 = 7.32 meses
