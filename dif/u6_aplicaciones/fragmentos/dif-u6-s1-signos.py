# ID: DIF-U6-S1
# Libro: Dif. U6 · subtema 6.1: puntos críticos e inflexión por cambio de signo
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-s1-signos.py
import sympy as sp

x = sp.symbols("x", real=True)
f = x**3 - 6 * x**2 + 9 * x + 1
d1, d2 = sp.diff(f, x), sp.diff(f, x, 2)
dx = sp.Rational(1, 10)            # menor que la distancia entre puntos
for c in sp.solve(d1, x):          # puntos críticos: f' = 0
    izq, der = d1.subs(x, c - dx), d1.subs(x, c + dx)
    tipo = ("máximo" if izq > 0 > der else
            "mínimo" if izq < 0 < der else "sin extremo")
    print(f"x = {c}: f' pasa de {izq} a {der} -> {tipo}")
for p in sp.solve(d2, x):          # candidatos a inflexión: f'' = 0
    izq, der = d2.subs(x, p - dx), d2.subs(x, p + dx)
    es = "inflexión" if izq * der < 0 else "sin inflexión"
    print(f"x = {p}: f'' pasa de {izq} a {der} -> {es}")
