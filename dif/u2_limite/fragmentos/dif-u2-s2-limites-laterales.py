# ID: DIF-U2-S2
# Libro: Dif. U2 · subtema 2.2: límites laterales leídos en una gráfica
# Repositorio: dif/u2_limite/fragmentos/dif-u2-s2-limites-laterales.py
import numpy as np
import matplotlib.pyplot as plt

def g(x):
    """Función de la figura: salto en -2, hueco en 1 y polo en 4."""
    if x < -2: return x + 3
    if x == 1: return 4.0
    if x < 4:  return 2 + (x - 1)**2 / 9
    if x > 4:  return 1 / (x - 4)
    return float("nan")                    # g(4) no está definida

for a in (-2, 1, 4):
    izq = [round(g(a - h), 4) for h in (0.1, 0.01, 0.001)]
    der = [round(g(a + h), 4) for h in (0.1, 0.01, 0.001)]
    print(f"a = {a:>2}: izquierda {izq}  derecha {der}  g(a) = {g(a)}")

xs = np.linspace(-4, 6, 2001)
ys = np.array([g(x) for x in xs])
ys[np.abs(ys) > 6] = np.nan                # no se dibuja cerca del polo
plt.plot(xs, ys, ".", ms=2)                # puntos sueltos, sin unir saltos
plt.ylim(-1.5, 6); plt.grid(True); plt.xlabel("x"); plt.show()
