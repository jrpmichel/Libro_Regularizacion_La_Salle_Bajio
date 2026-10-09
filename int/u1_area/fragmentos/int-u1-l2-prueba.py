# ID: INT-U1-L2
# Libro: Int. U1 · Laboratorio del Error 1.2 (prueba con código)
# Repositorio: int/u1_area/fragmentos/int-u1-l2-prueba.py
P = [3.0, 3.4, 3.6, 3.5, 3.1]       # kW, una lectura cada 10 s
dt = 10                             # s
print("sum(P)*dt =", round(sum(P) * dt, 6), "kJ  <- cinco rectángulos")
izq = sum(P[:-1]) * dt              # cuatro intervalos
der = sum(P[1:]) * dt
print("izquierda:", round(izq, 6), " derecha:", round(der, 6),
      " promedio:", round((izq + der) / 2, 6))
duracion = dt * (len(P) - 1)        # 40 s
print("cota gruesa:", round(max(P) * duracion, 6), "kJ")    # 144
