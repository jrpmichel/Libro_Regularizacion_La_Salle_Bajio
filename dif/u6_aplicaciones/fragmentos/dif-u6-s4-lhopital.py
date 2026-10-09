# ID: DIF-U6-S4
# Libro: Dif. U6 · subtema 6.4: L'Hopital dos veces en (1 - cos x)/x^2
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-s4-lhopital.py
import math
import sympy as sp

x = sp.symbols("x")
num, den = 1 - sp.cos(x), x**2
for paso in range(1, 3):
    cero = (num.subs(x, 0), den.subs(x, 0))
    print(f"paso {paso}: ({num})/({den}) -> en 0: {cero[0]}/{cero[1]}")
    num, den = sp.diff(num, x), sp.diff(den, x)
print("ya no es 0/0:", num.subs(x, 0), "/", den.subs(x, 0))   # 1/2
for h in (0.1, 0.01, 0.001):
    print(h, (1 - math.cos(h)) / h**2)                        # tiende a 0.5
