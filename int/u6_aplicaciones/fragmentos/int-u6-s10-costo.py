# ID: INT-U6-S10
# Libro: Int. U6 · subtema 6.10: costo desde el marginal y curva de aprendizaje
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-s10-costo.py
import sympy as sp

q = sp.symbols("q", positive=True)
Cm = 120 - sp.Rational(2, 5) * q             # costo marginal (pesos/pieza)
print("piezas 101 a 150:", sp.integrate(Cm, (q, 100, 150)), "pesos")
tq = 10 * q**sp.Rational(-32, 100)            # minutos por pieza
total = sp.integrate(tq, (q, 0, 100))
print("100 piezas:", round(float(total), 1), "min")         # 336.9
print("pieza 100:", round(float(tq.subs(q, 100)), 2), "min")  # 2.29
