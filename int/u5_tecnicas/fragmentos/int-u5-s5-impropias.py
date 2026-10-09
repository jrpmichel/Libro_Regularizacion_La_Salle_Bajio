# ID: INT-U5-S5
# Libro: Int. U5 · subtema 5.5: integrales parciales hacia infinito y hacia 0
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-s5-impropias.py
from scipy.integrate import quad

print("   R    1/x^2      1/x")
for R in (10, 100, 1000, 10000):
    a = quad(lambda x: 1 / x**2, 1, R, limit=200)[0]
    b = quad(lambda x: 1 / x, 1, R, limit=200)[0]
    print(f"{R:5d} {a:8.5f} {b:8.5f}")        # la primera se acerca a 1
print("  eps   1/sqrt(x)  1/x")
for e in (0.1, 0.01, 0.001, 0.0001):
    a = quad(lambda x: x**-0.5, e, 1)[0]
    b = quad(lambda x: 1 / x, e, 1)[0]
    print(f"{e:6} {a:8.5f} {b:8.5f}")         # la primera se acerca a 2
