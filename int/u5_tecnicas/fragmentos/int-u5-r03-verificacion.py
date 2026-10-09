# ID: INT-U5-R03
# Libro: Int. U5 · problema resuelto INT-U5-03 (verificación)
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-r03-verificacion.py
import sympy as sp

x = sp.symbols("x")
f = sp.sin(x) * sp.cos(x)
A = sp.sin(x)**2 / 2
respuestas = {"A": A, "B": -sp.cos(x)**2 / 2, "C": -sp.cos(2 * x) / 4}
for nombre, r in respuestas.items():
    app = sp.simplify(r - A) == 0                  # prueba de la aplicación
    buena = sp.simplify(sp.diff(r, x) - f) == 0    # prueba por la derivada
    print(nombre, " app:", app, " derivada:", buena,
          " r - A =", sp.simplify(r - A))
