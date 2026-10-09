# ID: DIF-U4-R04
# Libro: Dif. U4 · problema resuelto DIF-U4-04 (verificación)
# Repositorio: dif/u4_derivada/fragmentos/dif-u4-r04-verificacion.py
def iA(t):                         # rampa: i en A, t en s
    return 0.0 if t < 0 else (400 * t if t <= 0.005 else 2.0)


def iB(t):                         # escalón
    return 0.0 if t < 0 else 2.0


def laterales(i, t0, h=1e-6):
    izq = (i(t0) - i(t0 - h)) / h
    der = (i(t0 + h) - i(t0)) / h
    return round(izq, 3), round(der, 3)


print("A en 0:   ", laterales(iA, 0))          # (0.0, 400.0)
print("A en 5 ms:", laterales(iA, 0.005))      # (400.0, 0.0)
for h in (1e-3, 1e-6):
    print("B en 0:   ", laterales(iB, 0, h))   # izq = 2/h sin cota; der 0
print("v_A en la rampa:", 0.05 * 400, "V")     # L di/dt = 20 V
