# ID: DIF-U6-L1
# Libro: Dif. U6 · Laboratorio del Error 6.1 (prueba con código)
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-l1-prueba.py
import sympy as sp

x = sp.symbols("x", real=True)                  # cm
V = x * (60 - 2 * x) * (40 - 2 * x)             # cm^3
crit = [c for c in sp.solve(sp.diff(V, x), x) if 0 < c < 20]
print("crítico:", [round(float(c), 2) for c in crit])     # 7.85
print("V en el crítico:", round(float(V.subs(x, crit[0]))))   # 8450
for c in (0, 5):                                # la máquina: 0 <= x <= 5
    print(f"V({c}) = {V.subs(x, c)}")           # 0 y 7500
