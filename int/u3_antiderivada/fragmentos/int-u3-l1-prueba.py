# ID: INT-U3-L1
# Libro: Int. U3 · Laboratorio del Error 3.1 (prueba con código)
# Repositorio: int/u3_antiderivada/fragmentos/int-u3-l1-prueba.py
import sympy as sp

x = sp.symbols("x", positive=True)
f = x * sp.sqrt(x)
mal = x**sp.Rational(7, 2) / 3                # producto de las integrales
bien = sp.Rational(2, 5) * x**sp.Rational(5, 2)
print("derivada de la versión mal: ", sp.diff(mal, x))
print("derivada de la versión bien:", sp.diff(bien, x))
print("integrando:", sp.simplify(f))
print("en x = 4:", sp.diff(mal, x).subs(x, 4), "contra", f.subs(x, 4))
