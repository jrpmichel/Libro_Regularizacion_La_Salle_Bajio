# ID: DIF-U4-S2
# Libro: Dif. U4 · subtema 4.2: secantes del café con h cada vez menor
# Repositorio: dif/u4_derivada/fragmentos/dif-u4-s2-paso-al-limite.py
import math

T = lambda t: 20 + 70 * math.exp(-0.1 * t)    # temperatura del café, °C
t0 = 5                                         # minuto 5

print(" h (min)   pendiente de la secante (°C/min)")
for h in [5, 1, 0.1, 0.01, 0.001, 0.0001, -0.0001, -0.001, -0.01, -0.1,
          -1, -5]:
    m = (T(t0 + h) - T(t0)) / h
    print(f"{h:>8}   {m:.4f}")
# por la derecha y por la izquierda se acerca a -4.2457 °C/min
