# ID: DIF-U5-S5
# Libro: Dif. U5 · subtema 5.5: pendientes de e^x, 2^x y del seno (radianes y grados)
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-s5-exp-trig.py
import math

h = 1e-6
for a in (-1.0, 0.0, 1.0):
    m = (math.exp(a + h) - math.exp(a)) / h
    print(f"e^x en {a:4}: pendiente {m:.5f}   altura {math.exp(a):.5f}")
print("2^x en 0:", round((2**h - 1) / h, 5), "  ln 2 =", round(math.log(2), 5))

a = math.pi / 3                          # 60 grados, en radianes
m = (math.sin(a + h) - math.sin(a)) / h
print("seno en pi/3:", round(m, 5), "  cos =", round(math.cos(a), 5))

sen_g = lambda t: math.sin(math.radians(t))    # angulo t en grados
m = (sen_g(60 + h) - sen_g(60)) / h
print("seno en 60 grados, por grado:", round(m, 6),
      "  (pi/180)cos 60 =", round(math.pi / 180 * 0.5, 6))
