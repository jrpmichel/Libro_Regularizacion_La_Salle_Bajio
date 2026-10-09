# ID: DIF-U5-S2
# Libro: Dif. U5 · subtema 5.2: regla del producto contra el producto de derivadas
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-s2-producto.py
u, du = lambda x: x**2 + 1, lambda x: 2 * x
v, dv = lambda x: x**3 - x, lambda x: 3 * x**2 - 1
x0, h = 2.0, 1e-6

numerica = (u(x0 + h) * v(x0 + h) - u(x0) * v(x0)) / h
print("cociente de u*v :", round(numerica, 3))                  # 79.0
print("u'v + uv'       :", du(x0) * v(x0) + u(x0) * dv(x0))     # 79.0
print("u'v' (no sirve) :", du(x0) * dv(x0))                     # 44.0
