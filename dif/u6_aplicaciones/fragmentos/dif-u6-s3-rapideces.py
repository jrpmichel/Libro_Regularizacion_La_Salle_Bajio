# ID: DIF-U6-S3
# Libro: Dif. U6 · subtema 6.3: globo esférico que se infla a 100 cm^3/s
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-s3-rapideces.py
import sympy as sp

t = sp.symbols("t")
r = sp.Function("r")(t)                  # el radio depende del tiempo
V = sp.Rational(4, 3) * sp.pi * r**3
dVdt = sp.diff(V, t)                     # 4*pi*r^2 * dr/dt
drdt = sp.solve(sp.Eq(dVdt, 100), sp.diff(r, t))[0]
print("dr/dt =", drdt)
valor = drdt.subs(r, 5)                  # sustituir DESPUÉS de derivar
print("en r = 5 cm:", valor, "=", round(float(valor), 4), "cm/s")
