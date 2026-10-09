# ID: INT-U4-S3
# Libro: Int. U4 · subtema 4.3: integral con signo y área total
# Repositorio: int/u4_tfc/fragmentos/int-u4-s3-signo.py
import sympy as sp

x = sp.symbols("x", real=True)
f = sp.sin(x)
a, b = 0, 3 * sp.pi / 2
ceros = sorted(sp.solveset(f, x, sp.Interval.open(a, b)))     # [pi]
cortes = [a] + ceros + [b]
partes = [sp.integrate(f, (x, p, q)) for p, q in zip(cortes, cortes[1:])]
print("ceros:", ceros, "  integral de cada región:", partes)  # [2, -1]
print("integral con signo:", sum(partes))                     # 1
print("área total:", sum(abs(v) for v in partes))             # 3
