# ID: INT-U6-R10
# Libro: Int. U6 · problema resuelto INT-U6-10 (verificación)
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-r10-verificacion.py
import sympy as sp

t, T = sp.symbols("t T", nonnegative=True)
CA = 250000 + sp.integrate(20000 + 6000 * t, (t, 0, T))
CB = 400000 + sp.integrate(12000 + 1500 * t, (t, 0, T))
Tc = [s for s in sp.solve(sp.Eq(CA, CB), T) if s > 0][0]
print("costos iguales en T =", Tc, "=", round(float(Tc), 2), "años")
for Th in (5, 10):
    print(f"T = {Th}: A {CA.subs(T, Th)}, B {CB.subs(T, Th)} pesos")
