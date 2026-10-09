# ID: DIF-U2-L1
# Libro: Dif. U2 · Laboratorio del Error 2.1 (prueba con código)
# Repositorio: dif/u2_limite/fragmentos/dif-u2-l1-prueba.py
import sympy as sp

x = sp.symbols("x")
f = (x**2 - 9) / (x - 3)
print("sustitución directa:", f.subs(x, 3))    # nan: 0/0 no es un número
for h in (0.1, 0.01, 0.001):
    izq, der = f.subs(x, 3 - h), f.subs(x, 3 + h)
    print(3 - h, round(float(izq), 4), "|", 3 + h, round(float(der), 4))
print("límite:", sp.limit(f, x, 3))            # 6

# tres cocientes que empiezan como 0/0 en x = 0
for g in ((x**2 + 2*x) / x, (x**3 + x**2) / x, (x**2 + x) / x**2):
    print(g, "->", sp.limit(g, x, 0, "+"))     # 2, 0, oo
