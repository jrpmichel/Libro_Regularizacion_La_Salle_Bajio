# ID: DIF-U2-L2
# Libro: Dif. U2 · Laboratorio del Error 2.2 (prueba con código)
# Repositorio: dif/u2_limite/fragmentos/dif-u2-l2-prueba.py
import math

f = lambda x: (1 - math.cos(x)) / x**2            # límite exacto: 1/2
estable = lambda x: 2 * math.sin(x / 2)**2 / x**2  # misma función, sin resta
for k in range(1, 10):
    x = 10.0**-k
    print(f"x = 1e-{k}:  directa {f(x):.10f}   reescrita {estable(x):.10f}")
print(math.cos(1e-8) == 1.0)   # True: con 16 cifras, cos(1e-8) es 1 exacto
