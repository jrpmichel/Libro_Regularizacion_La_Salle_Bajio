# ID: INT-U5-R01
# Libro: Int. U5 · problema resuelto INT-U5-01 (verificación)
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-r01-verificacion.py
import sympy as sp

x = sp.symbols("x", positive=True)
sirven = {"a": (x**2 * sp.cos(x**3), sp.sin(x**3) / 3),
          "d": (sp.log(x) / x, sp.log(x)**2 / 2),
          "e": (sp.exp(sp.sqrt(x)) / sp.sqrt(x), 2 * sp.exp(sp.sqrt(x)))}
for k, (f, F) in sirven.items():
    print(k, "F' = f:", sp.simplify(sp.diff(F, x) - f) == 0)    # True
for k, f in {"b": sp.cos(x**3), "c": x * sp.cos(x**3)}.items():
    F = sp.integrate(f, x)
    especiales = {type(g).__name__ for g in F.atoms(sp.Function)}
    print(k, "sympy usa:", sorted(especiales))   # gamma, hyper...
