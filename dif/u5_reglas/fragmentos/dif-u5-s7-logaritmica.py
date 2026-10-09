# ID: DIF-U5-S7
# Libro: Dif. U5 · subtema 5.7: derivada de x^x por derivación logarítmica
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-s7-logaritmica.py
import math

f = lambda x: x**x
h = 1e-6
print("     x   cociente   x^x(ln x + 1)   x*x^(x-1)")
for a in (0.5, 1.0, 2.0, math.e):
    cociente = (f(a + h) - f(a)) / h
    log_dif = f(a) * (math.log(a) + 1)     # derivación logarítmica
    potencia = a * a**(a - 1)              # regla de la potencia, mal usada
    print(f"{a:6.4f} {cociente:10.4f} {log_dif:15.4f} {potencia:11.4f}")
