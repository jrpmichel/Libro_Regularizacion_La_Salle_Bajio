# ID: DIF-U4-S1
# Libro: Dif. U4 · subtema 4.1: pendiente de una recta y de una secante
# Repositorio: dif/u4_derivada/fragmentos/dif-u4-s1-pendientes.py
def pendiente(p, q):
    """Pendiente de la recta que pasa por p = (x1, y1) y q = (x2, y2)."""
    (x1, y1), (x2, y2) = p, q
    if x1 == x2:
        raise ValueError("x1 = x2: la recta es vertical")
    return (y2 - y1) / (x2 - x1)


recta = lambda x: 0.5 * x + 1
for a, b in [(1, 3), (4, 8), (-2, 10)]:      # en la recta: siempre 0.5
    print("recta  ", a, b, pendiente((a, recta(a)), (b, recta(b))))

f = lambda x: x**2
for b in [3, 2, 1.5, 1.1]:                   # secantes desde P = (1, 1)
    print("secante", 1, b, pendiente((1, f(1)), (b, f(b))))
