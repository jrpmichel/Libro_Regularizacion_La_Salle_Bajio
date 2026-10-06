# ID: DIF-U1-L1
# Libro: Dif. U1 · Laboratorio del Error 1.1 (prueba con código)
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-l1-prueba.py
f = lambda x: x**2 - 3*x
x, h = 4, 0.5
print("f(x+h)             =", f(x + h))                       # 6.75
print("f(x)+h             =", f(x) + h)                       # 4.5, lo que uso la escena
print("cociente verdadero =", (f(x + h) - f(x)) / h)          # 5.5
print("cociente escena    =", ((f(x) + h) - f(x)) / h)        # 1.0
