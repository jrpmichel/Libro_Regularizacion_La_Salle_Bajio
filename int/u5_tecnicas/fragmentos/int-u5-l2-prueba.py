# ID: INT-U5-L2
# Libro: Int. U5 · Laboratorio del Error 5.2 (prueba con código)
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-l2-prueba.py
import sympy as sp

x = sp.symbols("x")
mal = x * sp.sin(x) - sp.cos(x)
bien = x * sp.sin(x) + sp.cos(x)
print("derivada de la incorrecta:", sp.diff(mal, x))   # x*cos(x) + 2*sin(x)
print("derivada de la correcta:", sp.diff(bien, x))    # x*cos(x)
print("incorrecta en [0, pi]:", mal.subs(x, sp.pi) - mal.subs(x, 0))  # 2
print("integral:", sp.integrate(x * sp.cos(x), (x, 0, sp.pi)))        # -2
