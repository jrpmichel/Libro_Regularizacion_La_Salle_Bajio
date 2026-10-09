# ID: DIF-U5-L1
# Libro: Dif. U5 · Laboratorio del Error 5.1 (prueba con código)
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-l1-prueba.py
import math

w = 120 * math.pi                      # rad/s (60 Hz)
v = lambda t: 170 * math.sin(w * t)    # V
h = 1e-7
print("cociente en t = 0:", round((v(h) - v(0)) / h), "V/s")   # 64088
print("escena, 170cos(0):", 170, "V/s")       # falta la derivada interna
print("170*120*pi*cos(0):", round(170 * w), "V/s")   # 64088: cadena
