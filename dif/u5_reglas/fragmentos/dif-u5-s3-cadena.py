# ID: DIF-U5-S3
# Libro: Dif. U5 · subtema 5.3: la derivada de (x^2 + 1)^10 en x = 1 (semilla)
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-s3-cadena.py
f = lambda x: (x**2 + 1)**10
a = 1.0
for h in (1e-3, 1e-6):
    print(f"h = {h:g}: cociente = {(f(a + h) - f(a)) / h:.1f}")
# 5120: a la apuesta de la semilla le falta la derivada interna
print("10(x^2+1)^9      =", 10 * (a**2 + 1)**9)
# 10240: regla de la cadena
print("10(x^2+1)^9 * 2x =", 10 * (a**2 + 1)**9 * 2 * a)
