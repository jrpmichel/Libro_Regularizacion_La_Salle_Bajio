# ID: INT-U6-R05
# Libro: Int. U6 · problema resuelto INT-U6-05 (verificación)
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-r05-verificacion.py
import sympy as sp

x = sp.symbols("x", real=True)
lam = 4 - 4 * x                             # kg/m, x desde la articulación
L = sp.Rational(1, 2)
m = sp.integrate(lam, (x, 0, L))
xc = sp.integrate(x * lam, (x, 0, L)) / m
I0 = sp.integrate(x**2 * lam, (x, 0, L))
print("m =", m, "; x_c =", xc, "; I_O =", I0)   # 3/2, 2/9, 5/48
I = I0 + sp.Rational(1, 2) * L**2           # pinza de 0.5 kg en la punta
print("con pinza: I =", round(float(I), 4), "kg m^2; par =",
      round(float(6 * I), 3), "N m")         # 0.2292; 1.375
