# ID: INT-U6-S1
# Libro: Int. U6 · subtema 6.1: área entre sen x y cos x en [0, pi]
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-s1-area.py
import sympy as sp

x = sp.symbols("x", real=True)
f, g = sp.sin(x), sp.cos(x)
c = [r for r in sp.solve(sp.Eq(f, g), x) if 0 <= r <= sp.pi][0]
partes = [sp.integrate(f - g, (x, 0, c)),
          sp.integrate(f - g, (x, c, sp.pi))]
print("cruce en x =", c)                                 # pi/4
print("partes con signo:", partes)              # 1 - sqrt(2), 1 + sqrt(2)
print("integral con signo:", sp.simplify(sum(partes)))   # 2
print("área:", sp.simplify(sum(abs(p) for p in partes)))  # 2*sqrt(2)
