# ID: DIF-U2-S6
# Libro: Dif. U2 · subtema 2.6: clasificar una discontinuidad con sus laterales
# Repositorio: dif/u2_limite/fragmentos/dif-u2-s6-tipos-discontinuidad.py
import sympy as sp

x = sp.symbols("x")
f = (x**2 - 1) / (x**2 - 3*x + 2)       # = (x-1)(x+1) / ((x-1)(x-2))
for a in (1, 2):
    izq, der = sp.limit(f, x, a, "-"), sp.limit(f, x, a, "+")
    valor = f.subs(x, a)                # nan o zoo: f(a) no está definida
    if izq == der and izq.is_finite:
        tipo = "continua" if valor == izq else "evitable"
    elif izq.is_finite and der.is_finite:
        tipo = "de salto"
    else:
        tipo = "infinita"
    print(f"a = {a}: izq. {izq}, der. {der}, f(a) = {valor} -> {tipo}")
