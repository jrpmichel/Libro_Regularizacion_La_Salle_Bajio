# ID: DIF-U6-S2
# Libro: Dif. U6 · subtema 6.2: Newton contra bisección para x^2 - 2 = 0
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-s2-newton.py
f = lambda x: x**2 - 2
df = lambda x: 2 * x

tol = 1e-10
x, k, paso = 1.0, 0, 1.0               # Newton desde x0 = 1
while abs(f(x)) > tol or abs(paso) > tol:
    paso = f(x) / df(x)
    x, k = x - paso, k + 1
print(f"Newton:    {x:.12f} en {k} iteraciones")

a, b, k = 1.0, 2.0, 0                  # bisección en [1, 2] (Unidad 2)
while b - a > tol:
    m = (a + b) / 2
    a, b = (a, m) if f(a) * f(m) <= 0 else (m, b)
    k += 1
print(f"Bisección: {(a + b) / 2:.12f} en {k} iteraciones")
