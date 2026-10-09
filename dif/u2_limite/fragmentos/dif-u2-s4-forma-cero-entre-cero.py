# ID: DIF-U2-S4
# Libro: Dif. U2 · subtema 2.4: la forma 0/0 antes y después de simplificar
# Repositorio: dif/u2_limite/fragmentos/dif-u2-s4-forma-cero-entre-cero.py
import sympy as sp

x = sp.symbols("x")
f = (x**3 - 1) / (x - 1)
print(f.subs(x, 1))                  # nan: 0/0, todavía no hay respuesta
print(sp.factor(x**3 - 1))           # (x - 1)*(x**2 + x + 1)
g = sp.cancel(f)                     # cancela (x - 1): cerca de 1, x != 1
print(g, "->", g.subs(x, 1))         # x**2 + x + 1 -> 3
print(sp.limit(f, x, 1, "-"), sp.limit(f, x, 1, "+"))       # 3 3

r = (sp.sqrt(x + 1) - 1) / x         # con raíz: se usa el conjugado
print(sp.expand((sp.sqrt(x + 1) - 1) * (sp.sqrt(x + 1) + 1)))   # x
print(sp.limit(r, x, 0), "=", (1 / (sp.sqrt(x + 1) + 1)).subs(x, 0))
