# ID: DIF-U4-L2
# Libro: Dif. U4 · Laboratorio del Error 4.2 (prueba con código)
# Repositorio: dif/u4_derivada/fragmentos/dif-u4-l2-prueba.py
import sympy as sp

x = sp.symbols("x", real=True)
f, df, a = x**2, 2 * x, 3            # f'(x) = 2x (sección 4.3)
mal = f.subs(x, a) + df * (x - a)                 # pendiente f'(x)
bien = f.subs(x, a) + df.subs(x, a) * (x - a)     # pendiente f'(a)
print(sp.expand(mal), "|", sp.expand(bien))       # 2x²-6x+9 | 6x-9
xs = sp.Rational(31, 10)
print([float(e.subs(x, xs)) for e in (mal, bien, f)])  # 9.62 9.6 9.61
print(mal.subs(x, 4), bien.subs(x, 4))            # 17 15
print(sp.degree(mal, x), sp.degree(bien, x))      # 2 1: no es recta
