# ID: DIF-U3-R02
# Libro: Dif. U3 · problema resuelto DIF-U3-02 (verificación)
# Repositorio: dif/u3_transformaciones/fragmentos/dif-u3-r02-verificacion.py
import sympy as sp

x = sp.symbols("x", real=True)
a, b = sp.symbols("a b", positive=True)
f = lambda u: 4 - u**2                         # arco base: claro 4, flecha 4
sol = sp.solve([sp.Eq(b * 15, 2), sp.Eq(a * 4, 6)], [a, b])
g = sol[a] * f(sol[b] * x)
print(sol, "->", sp.expand(g))                 # a = 3/2, b = 2/15
print([g.subs(x, v) for v in (0, 5, 10, 15)])  # 6, 16/3, 10/3, 0
