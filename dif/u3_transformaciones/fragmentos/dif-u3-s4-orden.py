# ID: DIF-U3-S4
# Libro: Dif. U3 · subtema 3.4: factorizar b y el orden de las operaciones
# Repositorio: dif/u3_transformaciones/fragmentos/dif-u3-s4-orden.py
import sympy as sp

x = sp.symbols("x")
f = lambda u: u**2
b, c = 2, -6                                   # f(bx + c) = f(b(x - h))
h = sp.Rational(-c, b)
print("h =", h)                                # 3, no 6
print(sp.expand(f(b*x + c) - f(b*(x - h))))    # 0: es la misma función
print(sp.solve(sp.Eq(b*x + c, 0), x))          # vértice en [3]
print(sp.expand(f(b*(x - 6))))                 # desplazar 6 da otra función

s = sp.sin(x)
print(2*s + 3, "-> rango [1, 5]")
print(sp.expand(2*(s + 3)), "-> rango [4, 8]") # el orden vertical importa
