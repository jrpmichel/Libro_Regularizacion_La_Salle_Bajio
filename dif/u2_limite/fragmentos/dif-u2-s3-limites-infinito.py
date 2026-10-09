# ID: DIF-U2-S3
# Libro: Dif. U2 · subtema 2.3: límites en el infinito y valor final
# Repositorio: dif/u2_limite/fragmentos/dif-u2-s3-limites-infinito.py
import numpy as np
import sympy as sp

x = sp.symbols("x")
r = (3*x**2 + 2*x) / (x**2 + 4)
for v in (10, 100, 1000, 10**6):
    print(f"r({v}) = {float(r.subs(x, v)):.6f}   "
          f"r({-v}) = {float(r.subs(x, -v)):.6f}")
print("en +oo:", sp.limit(r, x, sp.oo), "| en -oo:", sp.limit(r, x, -sp.oo))
print("r(x) = 3 en x =", sp.solve(sp.Eq(r, 3), x))   # cruza su asíntota

tau = 0.5                                           # constante de tiempo (s)
for k in (1, 3, 5):
    print(f"t = {k} tau = {k*tau} s:  v = {12*(1 - np.exp(-k)):.3f} V de 12 V")
