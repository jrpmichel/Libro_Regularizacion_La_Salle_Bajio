# ID: DIF-U5-R01
# Libro: Dif. U5 · problema resuelto DIF-U5-01 (verificación)
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-r01-verificacion.py
import sympy as sp

x = sp.symbols("x", positive=True)
f = 4 * x**5 - 3 / x**2 + 7 * sp.sqrt(x) - 2
mia = 20 * x**4 + 6 / x**3 + sp.Rational(7, 2) / sp.sqrt(x)
print(sp.simplify(sp.diff(f, x) - mia) == 0)        # True
print(mia.subs(x, 1))                               # 59/2 = 29.5
fn = sp.lambdify(x, f)
print(round((fn(1 + 1e-6) - fn(1)) / 1e-6, 4))      # 29.5: cociente
