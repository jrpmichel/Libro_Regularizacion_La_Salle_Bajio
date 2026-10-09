# ID: INT-U3-S3
# Libro: Int. U3 · subtema 3.3: la condición y(2) = 5 escoge una curva
# Repositorio: int/u3_antiderivada/fragmentos/int-u3-s3-condicion.py
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

x, C = sp.symbols("x C")
f = 3 * x**2 - 3                         # y' = f(x)
x0, y0 = 2, 5                            # condición y(2) = 5
F = sp.integrate(f, x) + C
c = sp.solve(sp.Eq(F.subs(x, x0), y0), C)[0]
y = F.subs(C, c)
print("C =", c, "   y(x) =", y)          # C = 3

xx = np.linspace(-2.3, 2.6, 200)
for k in (-3, -1, 1, 5, 7):
    plt.plot(xx, sp.lambdify(x, F.subs(C, k))(xx), ":", color="gray")
plt.plot(xx, sp.lambdify(x, y)(xx), lw=2, label="y(2) = 5")
plt.plot(x0, y0, "o", color="red")
plt.xlabel("x"); plt.ylabel("y"); plt.legend(); plt.grid(True)
plt.show()
