# ID: INT-U2-S2
# Libro: Int. U2 · subtema 2.2: fórmulas cerradas de las sumas
# Repositorio: int/u2_riemann/fragmentos/int-u2-s2-sumas-cerradas.py
import sympy as sp

k, n = sp.symbols("k n", integer=True, positive=True)
for p in (1, 2, 3):
    cerrada = sp.factor(sp.summation(k**p, (k, 1, n)))   # fórmula cerrada
    directa = sum(j**p for j in range(1, 101))           # término a término
    print(f"suma de k^{p}: {cerrada};  n = 100: {cerrada.subs(n, 100)} = {directa}")
print("suma de 3, de k = 1 a n:", sp.summation(3, (k, 1, n)))   # 3n, no 3
