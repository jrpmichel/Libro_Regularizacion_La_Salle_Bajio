# ID: INT-U3-R01
# Libro: Int. U3 · problema resuelto INT-U3-01 (verificación)
# Repositorio: int/u3_antiderivada/fragmentos/int-u3-r01-verificacion.py
import sympy as sp

x = sp.symbols("x")
f = 6 * x**2 + 2
F = {"F1": 2*x**3 + 2*x, "F2": 2*x**3 + 2*x - 7, "F3": 6*x**3 + 2*x,
     "F4": 2*(x**3 + x + 1), "F5": 2*x**3 + 2*x + sp.pi}
for nombre, Fi in F.items():
    es = sp.simplify(sp.diff(Fi, x) - f) == 0
    resta = sp.simplify(Fi - F["F1"])
    print(nombre, "es antiderivada:", es, "   F - F1 =", resta)
