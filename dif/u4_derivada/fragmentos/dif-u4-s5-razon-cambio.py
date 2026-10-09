# ID: DIF-U4-S5
# Libro: Dif. U4 · subtema 4.5: velocidad media e instantánea
# Repositorio: dif/u4_derivada/fragmentos/dif-u4-s5-razon-cambio.py
s = lambda t: 4.9 * t**2        # caída libre: s en metros, t en segundos

for h in [1, 0.1, 0.01, 0.001]:
    v_media = (s(2 + h) - s(2)) / h
    print(f"[2, {2 + h}] s: velocidad media = {v_media:.4f} m/s")
# 4.9((2 + h)**2 - 4)/h = 4.9(4 + h)  ->  19.6 m/s en t = 2 s
