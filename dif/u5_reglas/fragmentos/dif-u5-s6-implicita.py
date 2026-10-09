# ID: DIF-U5-S6
# Libro: Dif. U5 · subtema 5.6: pendiente de la circunferencia x^2 + y^2 = 25
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-s6-implicita.py
import sympy as sp

x, y = sp.symbols("x y")
F = x**2 + y**2 - 25                # la curva es F(x, y) = 0
dydx = sp.idiff(F, y, x)            # derivación implícita
print("dy/dx =", dydx)              # -x/y
for px, py in [(3, 4), (3, -4), (0, 5)]:
    print(f"en ({px}, {py}): pendiente {dydx.subs({x: px, y: py})}")

rama = sp.sqrt(25 - x**2)           # rama superior, despejada
print("rama superior en x = 3:", sp.diff(rama, x).subs(x, 3))   # -3/4
