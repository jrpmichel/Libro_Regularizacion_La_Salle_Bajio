# ID: DIF-U6-L3
# Libro: Dif. U6 · Laboratorio del Error 6.3 (prueba con código)
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-l3-prueba.py
import math

dV = -0.4                                   # m^3/min
h = 2.0                                     # m; entonces r = h/2 = 1 m
escena = 3 * dV / (math.pi * 1.0**2)        # r fijo en 1 m antes de derivar
bien = 4 * dV / (math.pi * h**2)            # V = pi h^3/12, derivado
print("escena  :", round(escena, 4), "m/min")    # -0.382
print("correcta:", round(bien, 4), "m/min")      # -0.127
dt = 1e-4                                   # comprobación: nivel tras dt
h2 = (h**3 + 12 * dV * dt / math.pi)**(1 / 3)   # V(h2) = V(h) + dV*dt
print("numérica:", round((h2 - h) / dt, 4), "m/min")
