# ID: INT-U4-S4
# Libro: Int. U4 · subtema 4.4: valor medio y los puntos c
# Repositorio: int/u4_tfc/fragmentos/int-u4-s4-valor-medio.py
import numpy as np
import sympy as sp
from scipy.optimize import brentq

x = sp.symbols("x", real=True)
f = 2 + sp.sin(x)
a, b = 0, 4
fbar = sp.integrate(f, (x, a, b)) / (b - a)
print("valor medio:", sp.simplify(fbar), "=", round(float(fbar), 4))

g = sp.lambdify(x, f - fbar)                 # se busca f(c) = valor medio
xs = np.linspace(a, b, 401)
cs = [brentq(g, p, q) for p, q in zip(xs, xs[1:]) if g(p) * g(q) < 0]
print("puntos c:", [round(c, 4) for c in cs])   # 0.4262 y 2.7154
