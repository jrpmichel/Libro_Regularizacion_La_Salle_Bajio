# ID: INT-U1-L3
# Libro: Int. U1 · Laboratorio del Error 1.3 (prueba con código)
# Repositorio: int/u1_area/fragmentos/int-u1-l3-prueba.py
v_max = 90                          # km/h
T_min = 4                           # min: 2 acelerando y 2 frenando
print("con 4 min = 0.4 h:", 0.5 * 0.4 * v_max, "km")        # 18: mal
T_h = T_min / 60                    # 4 min = 0.0667 h
print("distancia:", round(0.5 * T_h * v_max, 6), "km")      # 3 km
print("cota (v_max * T):", round(v_max * T_h, 6), "km")     # 6 km
print("18 km en 4 min serían", round(18 / T_h, 6), "km/h")  # 270 > 90
