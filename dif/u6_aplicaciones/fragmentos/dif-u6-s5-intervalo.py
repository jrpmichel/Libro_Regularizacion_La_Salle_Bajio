# ID: DIF-U6-S5
# Libro: Dif. U6 · subtema 6.5: extremos absolutos en un intervalo cerrado
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-s5-intervalo.py
import sympy as sp

x = sp.symbols("x", real=True)
f = x**3 - 6 * x**2 + 9 * x + 1
a, b = sp.Rational(1, 2), sp.Rational(9, 2)
criticos = [c for c in sp.solve(sp.diff(f, x), x) if a <= c <= b]
candidatos = [a, b] + criticos               # extremos y puntos críticos
for c in candidatos:
    print(f"f({c}) = {f.subs(x, c)}")
valores = {c: f.subs(x, c) for c in candidatos}
print("máximo en x =", max(valores, key=valores.get))     # 9/2, extremo
print("mínimo en x =", min(valores, key=valores.get))     # 3, interior
