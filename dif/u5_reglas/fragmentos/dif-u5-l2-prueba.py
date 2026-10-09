# ID: DIF-U5-L2
# Libro: Dif. U5 · Laboratorio del Error 5.2 (prueba con código)
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-l2-prueba.py
import sympy as sp

x, y = sp.symbols("x y")
F = x**2 + x * y + y**2 - 7
print("correcta:", sp.idiff(F, y, x).subs({x: 1, y: 2}))     # -4/5
print("escena  :", sp.Rational(-2 * 1, 1 + 2 * 2))           # -2/5
rama = (-x + sp.sqrt(28 - 3 * x**2)) / 2      # y despejada; pasa por (1, 2)
m = sp.simplify(sp.diff(rama, x).subs(x, 1))
print("rama en x = 1:", rama.subs(x, 1), "| pendiente:", m)  # 2 | -4/5
