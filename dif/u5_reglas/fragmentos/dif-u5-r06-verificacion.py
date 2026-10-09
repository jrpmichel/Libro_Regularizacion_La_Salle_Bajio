# ID: DIF-U5-R06
# Libro: Dif. U5 · problema resuelto DIF-U5-06 (verificación)
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-r06-verificacion.py
import sympy as sp

V, P = sp.symbols("V P", positive=True)
C2 = 100 * 2**sp.Rational(7, 5)
iso = sp.idiff(P * V - 200, P, V)                  # P*V = 100 kPa * 2 L
adi = sp.idiff(P * V**sp.Rational(7, 5) - C2, P, V)
print(iso, "|", sp.simplify(adi))                  # -P/V | -7*P/(5*V)
punto = {P: 100, V: 2}
print(iso.subs(punto), adi.subs(punto))            # -50 y -70 kPa/L
print(100 + (-50) * (-0.1), 100 + (-70) * (-0.1))  # 105.0 y 107.0
print(round(200 / 1.9, 2), round(100 * (2 / 1.9)**1.4, 2))  # 105.26 107.45
