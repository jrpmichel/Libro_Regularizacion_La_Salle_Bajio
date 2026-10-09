# ID: INT-U5-L3
# Libro: Int. U5 · Laboratorio del Error 5.3 (prueba con código)
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-l3-prueba.py
import sympy as sp

x = sp.symbols("x")
e = sp.symbols("epsilon", positive=True)
derecha = sp.integrate(1 / x, (x, e, 1))                 # -log(epsilon)
print("parte derecha:", derecha, "-> límite:", sp.limit(derecha, e, 0, "+"))
izquierda = sp.integrate(1 / x, (x, -1, -e))
print("ventana (-e, e):", sp.simplify(izquierda + derecha))       # 0
otra = izquierda + sp.integrate(1 / x, (x, 2 * e, 1))
print("ventana (-e, 2e):", sp.simplify(otra))                     # -log(2)
