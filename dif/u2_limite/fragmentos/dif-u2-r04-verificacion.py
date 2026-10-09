# ID: DIF-U2-R04
# Libro: Dif. U2 · problema resuelto DIF-U2-04 (verificación)
# Repositorio: dif/u2_limite/fragmentos/dif-u2-r04-verificacion.py
def P(e):
    """Potencia del calentador (% de la nominal) según el error e (°C)."""
    if e <= 0: return 0.0
    if e < 5:  return 10 + 18 * e
    return 100.0

for a in (0, 5):
    izq = [round(P(a - h), 3) for h in (0.1, 0.01, 0.001)]
    der = [round(P(a + h), 3) for h in (0.1, 0.01, 0.001)]
    print(f"e = {a}: izquierda {izq}  derecha {der}  P({a}) = {P(a)}")
print("e de -0.1 a 0.1 °C:", P(-0.1), "->", round(P(0.1), 1), "%")
print("ley continua 20e:  ", 0.0, "->", 20 * 0.1, "%")
