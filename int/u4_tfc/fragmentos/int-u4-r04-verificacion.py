# ID: INT-U4-R04
# Libro: Int. U4 · problema resuelto INT-U4-04 (verificación)
# Repositorio: int/u4_tfc/fragmentos/int-u4-r04-verificacion.py
import numpy as np
import sympy as sp

x, eps = sp.symbols("x eps", positive=True)
F = -1 / x
print("F(1) - F(-1) =", F.subs(x, 1) - F.subs(x, -1))   # -2, sin sentido
f = lambda v: 1 / v**2
for n in (10, 100, 1000):                     # n par: ningún punto en 0
    xs = np.linspace(-1, 1, n + 1)
    M = f((xs[:-1] + xs[1:]) / 2).sum() * 2 / n
    print(f"M_{n} = {M:.1f}")
lado = sp.integrate(1 / x**2, (x, eps, 1))    # 1/eps - 1
print("de eps a 1:", lado)
print("con eps = 0.01:", lado.subs(eps, sp.Rational(1, 100)))   # 99
