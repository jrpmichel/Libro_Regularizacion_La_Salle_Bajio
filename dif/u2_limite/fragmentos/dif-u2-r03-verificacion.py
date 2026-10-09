# ID: DIF-U2-R03
# Libro: Dif. U2 · problema resuelto DIF-U2-03 (verificación)
# Repositorio: dif/u2_limite/fragmentos/dif-u2-r03-verificacion.py
import sympy as sp

x = sp.symbols("x")
f = (sp.sqrt(x + 3) - 2) / (x**2 - 1)
g = 1 / ((x + 1) * (sp.sqrt(x + 3) + 2))       # forma simplificada, x != 1
print("f(1):", f.subs(x, 1))                   # nan: 0/0
print("f - g:", sp.simplify(f - g))            # 0
print("laterales:", sp.limit(f, x, 1, "-"), sp.limit(f, x, 1, "+"))
print("g(1) =", g.subs(x, 1))                  # 1/8
for v in (0.99, 0.999, 1.001, 1.01):
    print(v, round(float(f.subs(x, v)), 5))    # cerca de 0.125
