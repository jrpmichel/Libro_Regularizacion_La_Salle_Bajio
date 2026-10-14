# ID: PRE-U0-S3
# Libro: Anexo C · sympy: derivadas, integrales, limites, ecuaciones y matrices
# Repositorio: pre/u0_primer_notebook/fragmentos/pre-u0-s3-sympy.py
import sympy as sp

x, h = sp.symbols("x h")
f = x**2 * sp.sin(x)
print("f'(x) =", sp.diff(f, x))
I = sp.integrate(sp.exp(-x), (x, 0, sp.oo))         # integral impropia
print("integral de e^(-x) en [0, oo):", I)
print("limite =", sp.limit(((x + h)**2 - x**2) / h, h, 0))
print("x**2 = 9:", sp.solve(x**2 - 9, x))          # las dos raices
print(sp.Rational(1, 10) + sp.Rational(2, 10), "contra", 0.1 + 0.2)
M = sp.Matrix([[2, 1], [1, 3]])
print("M^-1 =", M.inv(), " det =", M.det(), " rango =", M.rank())
g = sp.lambdify(x, f, "numpy")                     # f como funcion de numpy
print("f(pi/2) =", f.subs(x, sp.pi / 2), "  f(1.0) =", round(g(1.0), 6))
