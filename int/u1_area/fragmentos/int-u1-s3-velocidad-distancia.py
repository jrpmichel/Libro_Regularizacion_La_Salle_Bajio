# ID: INT-U1-S3
# Libro: Int. U1 · subtema 1.3: desplazamiento y distancia desde una tabla de velocidad
# Repositorio: int/u1_area/fragmentos/int-u1-s3-velocidad-distancia.py
t = [0.0, 4.0]                      # s
v = [4.0, -4.0]                     # m/s: v = 4 - 2t (figura 1.4b)
desp, dist = 0.0, 0.0
for k in range(len(t) - 1):
    t0, t1, v0, v1 = t[k], t[k + 1], v[k], v[k + 1]
    if v0 * v1 < 0:                 # cambia de signo: buscar el cruce
        tc = t0 + v0 * (t1 - t0) / (v0 - v1)
        partes = [(t0, tc, v0, 0.0), (tc, t1, 0.0, v1)]
    else:
        partes = [(t0, t1, v0, v1)]
    for a, b, va, vb in partes:
        area = (va + vb) / 2 * (b - a)   # trapecio con signo
        desp += area
        dist += abs(area)
print("desplazamiento:", desp, "m")      # 0.0
print("distancia recorrida:", dist, "m") # 8.0
