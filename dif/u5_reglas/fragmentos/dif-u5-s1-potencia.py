# ID: DIF-U5-S1
# Libro: Dif. U5 · subtema 5.1: regla de la potencia contra el cociente incremental
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-s1-potencia.py
# Regla de la potencia contra el cociente incremental en a = 3
a, h = 3.0, 1e-6
for n in [2, 3, 5, -1, -2, 0.5, 1.5]:
    cociente = ((a + h)**n - a**n) / h
    regla = n * a**(n - 1)
    print(f"n = {n:>4}: cociente {cociente:10.5f}   n*a^(n-1) {regla:10.5f}")
