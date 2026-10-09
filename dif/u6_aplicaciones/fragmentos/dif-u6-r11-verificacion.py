# ID: DIF-U6-R11
# Libro: Dif. U6 · problema resuelto DIF-U6-11 (verificación)
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-r11-verificacion.py
import sympy as sp

x = sp.symbols("x", real=True)                 # m
M = 12 * x - x**3 / 3                          # kN*m
V = sp.diff(M, x)                              # kN
print("V =", V, "| dV/dx =", sp.diff(V, x), "= -w(x)")
print("M(0), M(6):", M.subs(x, 0), M.subs(x, 6), "| V(0), V(6):",
      V.subs(x, 0), V.subs(x, 6))
xm = [c for c in sp.solve(V, x) if 0 <= c <= 6][0]
print("V = 0 en x =", xm, "=", round(float(xm), 3), "m")
print("M máx =", round(float(M.subs(x, xm)), 2), "| M(3) =", M.subs(x, 3),
      "| M(4) =", round(float(M.subs(x, 4)), 2))
