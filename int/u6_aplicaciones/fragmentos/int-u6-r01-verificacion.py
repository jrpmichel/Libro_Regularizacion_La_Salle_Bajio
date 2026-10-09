# ID: INT-U6-R01
# Libro: Int. U6 · problema resuelto INT-U6-01 (verificación)
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-r01-verificacion.py
import sympy as sp

t = sp.symbols("t", real=True)
P = 180 + 45 * t - 5 * t**2                     # kW, t en h
t1, t2 = sorted(sp.solve(sp.Eq(P, 250), t))     # 2 y 7
E = sp.integrate(P - 250, (t, t1, t2))
print("de", t1, "a", t2, "h: E =", E, "=", round(float(E), 2), "kWh")
tp = sp.solve(sp.diff(P, t), t)[0]
print("pico:", P.subs(t, tp), "kW en t =", tp, "h")      # 1125/4; 9/2
print("capacidad con 80 %:", round(float(E) / 0.8, 1), "kWh")  # 130.2
print("con signo en [0, 8]:", sp.integrate(P - 250, (t, 0, 8)))  # 80/3
