# ID: DIF-U6-R04
# Libro: Dif. U6 · problema resuelto DIF-U6-04 (verificación)
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-r04-verificacion.py
import math
import sympy as sp

n, s = sp.symbols("n s", positive=True)
cociente = 50 * n * sp.log(n) / n**sp.Rational(6, 5)
print("límite:", sp.limit(cociente, n, sp.oo))         # 0
print("L'Hopital:", sp.simplify(sp.diff(50 * sp.log(n), n)
                                / sp.diff(n**sp.Rational(1, 5), n)))
cruce = sp.nsolve(50 * s - sp.exp(s / 5), s, 37)       # n = e^s
print("cruce en n =", f"{math.exp(float(cruce)):.2e}")  # 2.38e+16
for N in (1e6, 1e9):
    print(f"n = {N:.0e}: T1/T2 = {50 * math.log(N) / N**0.2:.1f}")
