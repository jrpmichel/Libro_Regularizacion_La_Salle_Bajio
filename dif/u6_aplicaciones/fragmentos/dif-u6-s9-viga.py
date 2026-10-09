# ID: DIF-U6-S9
# Libro: Dif. U6 · subtema 6.9: viga simplemente apoyada con carga uniforme
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-s9-viga.py
import sympy as sp

x = sp.symbols("x", real=True)
w, L = 10, 8                                # kN/m, m
M = w * L * x / 2 - w * x**2 / 2            # momento, kN*m
V = sp.diff(M, x)                           # dM/dx = V
print("V(x) =", V, "| dV/dx =", sp.diff(V, x), "= -w")
xm = sp.solve(V, x)[0]                      # el máximo de M está donde V = 0
print("V = 0 en x =", xm, "m; M máximo =", M.subs(x, xm), "kN*m")
print("wL^2/8 =", w * L**2 / 8)
