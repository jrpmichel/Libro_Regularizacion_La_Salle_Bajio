# ID: INT-U6-L3
# Libro: Int. U6 · Laboratorio del Error 6.3 (prueba con código)
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-l3-prueba.py
# Sección T: patín 12 x 2 cm (centro en y = 11), alma 2 x 10 cm (y = 5)
partes = [(12, 2, 11), (2, 10, 5)]          # (b, h, y_c) en cm
A = sum(b * h for b, h, _ in partes)
yc = sum(b * h * y for b, h, y in partes) / A
I_ia = sum(b * h**3 / 12 for b, h, _ in partes)
I = sum(b * h**3 / 12 + b * h * (y - yc)**2 for b, h, y in partes)
print(f"centroide a {yc:.3f} cm de la base")         # 8.273
print(f"sin ejes paralelos: {I_ia:.1f} cm^4")         # 174.7
print(f"con ejes paralelos: {I:.1f} cm^4")            # 567.4
