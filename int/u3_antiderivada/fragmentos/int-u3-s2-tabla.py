# ID: INT-U3-S2
# Libro: Int. U3 · subtema 3.2: comprobación de la tabla de integrales
# Repositorio: int/u3_antiderivada/fragmentos/int-u3-s2-tabla.py
import numpy as np
import sympy as sp

x = sp.symbols("x", positive=True)       # rama x > 0
tabla = [(x**3, x**4 / 4), (1 / x, sp.log(x)), (sp.exp(x), sp.exp(x)),
         (2**x, 2**x / sp.log(2)), (sp.sin(x), -sp.cos(x)),
         (sp.cos(x), sp.sin(x)), (1 / sp.cos(x)**2, sp.tan(x)),
         (sp.cos(3 * x), sp.sin(3 * x) / 3)]      # argumento lineal
for f, F in tabla:
    ok = sp.simplify(sp.diff(F, x) - f) == 0
    print(f"{str(f):12s} -> {str(F):14s} derivada correcta: {ok}")

# rama x < 0: pendiente de ln|x| en x = -2 con una diferencia centrada
h = 1e-6
G = lambda t: np.log(np.abs(t))
print("pendiente en x = -2:", (G(-2 + h) - G(-2 - h)) / (2 * h))  # -0.5
