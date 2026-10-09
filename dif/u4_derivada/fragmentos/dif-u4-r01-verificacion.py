# ID: DIF-U4-R01
# Libro: Dif. U4 · problema resuelto DIF-U4-01 (verificación)
# Repositorio: dif/u4_derivada/fragmentos/dif-u4-r01-verificacion.py
C = lambda q: 2000 + 15 * q + 0.05 * q**2    # pesos, q piezas por día

m_sec = (C(120) - C(100)) / (120 - 100)      # secante de 100 a 120
m_org = C(100) / 100                         # recta desde el origen
print(C(100), C(120))                        # 4000.0 4520.0
print("secante:", m_sec, "pesos/pieza")      # 26.0
print("promedio de las 100:", m_org)         # 40.0
print([round(C(q + 1) - C(q), 2) for q in (100, 119)])  # 25.05 26.95
