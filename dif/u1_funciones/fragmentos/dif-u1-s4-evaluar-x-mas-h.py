# ID: DIF-U1-S4
# Libro: Dif. U1 · subtema 1.4: f(x+h) frente a f(x)+h
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-s4-evaluar-x-mas-h.py
import sympy as sp

x, h = sp.symbols("x h")
f = lambda u: u**2 - 3*u
print(sp.expand(f(x + h)))                  # f(x+h): sustituye x por (x+h) en todas partes
print(sp.simplify((f(x + h) - f(x)) / h))   # cociente de diferencias: h + 2*x - 3
print(sp.simplify(f(x + h) - (f(x) + h)))   # h*(h + 2*x - 4): no es 0, f(x+h) != f(x)+h
