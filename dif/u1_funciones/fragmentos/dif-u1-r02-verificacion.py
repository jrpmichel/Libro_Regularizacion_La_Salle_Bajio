# ID: DIF-U1-R02
# Libro: Dif. U1 · problema resuelto DIF-U1-02 (verificación)
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-r02-verificacion.py
import sympy as sp

x, y = sp.symbols("x y", real=True)
print(sp.solve(sp.Eq(3**2 + y**2, 25), y))        # [-4, 4]: dos puntos en x = 3
techo = sp.sqrt(25 - x**2)
print(techo.subs(x, 3), techo.subs(x, 0), techo.subs(x, 5))   # 4 5 0
