# ID: DIF-U3-L1
# Libro: Dif. U3 · Laboratorio del Error 3.1 (prueba con código)
# Repositorio: dif/u3_transformaciones/fragmentos/dif-u3-l1-prueba.py
f = lambda x: x**2
g = lambda x: f(x + 3)                         # (x + 3)^2

print("g(3) =", g(3), "| g(-3) =", g(-3))    # 36 | 0: vértice en -3
for x in (-5, -4, -3, -2, -1):
    print(x, g(x))                             # simétrica respecto de x = -3
