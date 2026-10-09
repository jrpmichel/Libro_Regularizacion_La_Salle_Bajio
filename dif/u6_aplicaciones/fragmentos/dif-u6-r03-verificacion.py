# ID: DIF-U6-R03
# Libro: Dif. U6 · problema resuelto DIF-U6-03 (verificación)
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-r03-verificacion.py
import sympy as sp

t = sp.symbols("t")
h = sp.Function("h")(t)                      # nivel, m
V = sp.pi * (h / 2)**2 * h / 3               # cono con r = h/2
dVdt = sp.Rational(-2, 5)                    # -0.4 m^3/min
dhdt = sp.solve(sp.Eq(sp.diff(V, t), dVdt), sp.diff(h, t))[0]
print("dh/dt =", sp.simplify(dhdt))          # -8/(5 pi h^2) = -1.6/(pi h^2)
for nivel in (2, 1):
    valor = float(dhdt.subs(h, nivel))
    print(f"h = {nivel} m: dh/dt = {valor:.4f} m/min")
