# ID: DIF-U2-L4
# Libro: Dif. U2 · Laboratorio del Error 2.4 (prueba con código)
# Repositorio: dif/u2_limite/fragmentos/dif-u2-l4-prueba.py
import sympy as sp
from sympy.calculus.util import continuous_domain

def biseccion(f, a, b, tol=1e-6):    # la de la escena: no revisa f(c)
    pasos = 0
    while b - a > tol:
        c = (a + b) / 2
        if f(a) * f(c) <= 0: b = c
        else: a = c
        pasos += 1
    return (a + b) / 2, pasos

f = lambda x: x / (x - 2) - 3
c, n = biseccion(f, 0, 2.5)
print(f"'raíz' aprox. {c:.6f} en {n} pasos, pero f(c) = {f(c):.3e}")
x = sp.symbols("x")
I = sp.Interval(0, sp.Rational(5, 2))
print(continuous_domain(x / (x - 2) - 3, x, I))     # falta el 2
print(sp.solve(sp.Eq(x / (x - 2) - 3, 0), x))       # [3]
c, n = biseccion(f, 2.5, 4)
print(f"en [2.5, 4]: raíz aprox. {c:.6f}, f(c) = {f(c):.1e}")
