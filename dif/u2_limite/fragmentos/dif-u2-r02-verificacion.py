# ID: DIF-U2-R02
# Libro: Dif. U2 · problema resuelto DIF-U2-02 (verificación)
# Repositorio: dif/u2_limite/fragmentos/dif-u2-r02-verificacion.py
P, L, a = 12.0, 6.0, 2.0                  # kN, m, m
RA, RB = P * (L - a) / L, P * a / L       # reacciones por equilibrio
V = lambda x: RA if x < a else RA - P     # cortante, 0 < x < L, x != a

izq = [V(a - h) for h in (0.1, 0.01, 0.001)]
der = [V(a + h) for h in (0.1, 0.01, 0.001)]
print("RA =", RA, "kN, RB =", RB, "kN, suma =", RA + RB)   # 8.0 4.0 12.0
print("izquierda:", izq, "| derecha:", der)                # 8.0 | -4.0
print("salto:", izq[-1] - der[-1], "kN")                   # 12.0 = P
