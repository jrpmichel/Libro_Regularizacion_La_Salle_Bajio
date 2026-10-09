# ID: DIF-U2-S7
# Libro: Dif. U2 · subtema 2.6: método de bisección
# Repositorio: dif/u2_limite/fragmentos/dif-u2-s7-biseccion.py
def biseccion(f, a, b, tol=1e-6):
    """Raíz de f en [a, b]. Exige f continua en [a, b] y f(a)·f(b) < 0."""
    if f(a) * f(b) >= 0:
        raise ValueError("f(a) y f(b) deben tener signos opuestos")
    pasos = 0
    while b - a > tol:
        c = (a + b) / 2
        if f(a) * f(c) <= 0: b = c           # la raíz queda en [a, c]
        else: a = c                          # la raíz queda en [c, b]
        pasos += 1
    c = (a + b) / 2
    return c, pasos

f = lambda x: x**3 - x - 2
c, n = biseccion(f, 1, 2)
print(f"raíz aprox. {c:.6f} en {n} pasos;  f(c) = {f(c):.1e}")
