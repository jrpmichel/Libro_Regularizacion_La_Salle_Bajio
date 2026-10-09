# ID: DIF-U6-R02
# Libro: Dif. U6 · problema resuelto DIF-U6-02 (verificación)
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-r02-verificacion.py
import sympy as sp

d = sp.Integer(36)                          # km entre casetas
for minutos in (18, sp.Rational(103, 5)):   # 10:00 a 10:18 y 20.6 min
    media = d / (minutos / 60)              # km/h
    print(f"{float(minutos)} min: media {float(media):.1f} km/h,",
          "mayor que 110:", media > 110)
# 120.0, True: hay un instante a 120 | 104.9, False: no concluye
