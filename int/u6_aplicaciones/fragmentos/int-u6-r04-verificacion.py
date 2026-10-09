# ID: INT-U6-R04
# Libro: Int. U6 · problema resuelto INT-U6-04 (verificación)
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-r04-verificacion.py
import math
import sympy as sp

y, b = sp.symbols("y b_f", positive=True)

def inercia(bf, tf, tw=sp.Rational(3, 5), h=30):
    """I centroidal (cm^4) de una sección I: integral de y^2 w(y) dy."""
    a = sp.Rational(h, 2) - tf
    return (2 * sp.integrate(y**2 * bf, (y, a, sp.Rational(h, 2)))
            + sp.integrate(y**2 * tw, (y, -a, a)))

for tf in (1, sp.Rational(13, 10)):
    bmin = sp.solve(sp.Eq(inercia(b, tf), 9000), b)[0]
    bf = math.ceil(float(bmin))                  # ancho entero en cm
    A = 2 * bf * tf + sp.Rational(3, 5) * (30 - 2 * tf)
    print(f"t_f = {float(tf)}: b_f >= {float(bmin):.2f}, con {bf} cm"
          f" I = {float(inercia(bf, tf)):.1f}, A = {float(A):.2f} cm^2")
print("rectángulo 4 x 30: I =", 4 * 30**3 / 12, ", A =", 4 * 30)
