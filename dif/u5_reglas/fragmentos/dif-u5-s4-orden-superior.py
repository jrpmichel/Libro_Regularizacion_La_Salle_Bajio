# ID: DIF-U5-S4
# Libro: Dif. U5 · subtema 5.4: derivadas sucesivas de un polinomio
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-s4-orden-superior.py
import sympy as sp

x = sp.symbols("x")
f = 2 * x**5 - x**3 + 4 * x
for n in range(1, 7):
    f = sp.diff(f, x)                 # cada vuelta deriva la derivada anterior
    print(f"derivada {n}: {f}")
