# ID: INT-U3-S1
# Libro: Int. U3 · subtema 3.1: familia de antiderivadas de f(x) = 2x
# Repositorio: int/u3_antiderivada/fragmentos/int-u3-s1-familia.py
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp

x = sp.symbols("x")
f = 2 * x
F = sp.integrate(f, x)                  # sympy da una antiderivada, sin C
print("F(x) =", F, "+ C")
print("F'(x) - f(x) =", sp.simplify(sp.diff(F, x) - f))    # 0

Fn = sp.lambdify(x, F, "numpy")
xx = np.linspace(-2, 2, 200)
for c in (-2, -1, 0, 1, 2):
    plt.plot(xx, Fn(xx) + c, label=f"C = {c}")
    plt.plot(1, Fn(1) + c, "o", color="red")   # pendiente 2 en x = 1
plt.xlabel("x"); plt.ylabel("x**2 + C")
plt.legend(); plt.grid(True); plt.show()
