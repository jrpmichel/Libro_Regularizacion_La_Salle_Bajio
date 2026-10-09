# ID: INT-U5-S3
# Libro: Int. U5 · subtema 5.3: descomposición de 2/(x^2 - 1)
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-s3-fracciones.py
import sympy as sp

x = sp.symbols("x")
f = 2 / (x**2 - 1)
piezas = sp.apart(f)                     # -1/(x + 1) + 1/(x - 1)
print("piezas:", piezas)
print("suman f:", sp.simplify(piezas - f) == 0)              # True
F = sp.integrate(piezas, x)
print("F =", F)                          # log(x - 1) - log(x + 1)
# sympy omite |·|: escribe ln(x - 1), válido para x > 1, como en [2, 3]
valor = F.subs(x, 3) - F.subs(x, 2)
print("de 2 a 3:", sp.simplify(valor), "=", round(float(valor), 6))
