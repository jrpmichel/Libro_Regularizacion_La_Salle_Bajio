# ID: DIF-U4-L3
# Libro: Dif. U4 · Laboratorio del Error 4.3 (prueba con código)
# Repositorio: dif/u4_derivada/fragmentos/dif-u4-l3-prueba.py
f = abs

for h in [0.1, 0.01, 0.001]:
    adelante = (f(0 + h) - f(0)) / h
    atras = (f(0) - f(0 - h)) / h
    centrada = (f(0 + h) - f(0 - h)) / (2 * h)
    print(h, adelante, atras, centrada)      # 1.0 -1.0 0.0 con cada h
# la centrada promedia 1 y -1; como esas dos no se acercan entre sí,
# f'(0) no existe aunque el programa imprima 0
