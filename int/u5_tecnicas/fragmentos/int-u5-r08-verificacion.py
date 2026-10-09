# ID: INT-U5-R08
# Libro: Int. U5 · problema resuelto INT-U5-08 (verificación)
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-r08-verificacion.py
import sympy as sp

t, R = sp.symbols("t R", positive=True)
f = 2 / sp.pi / (1 + t**2)
print("integral de f:", sp.integrate(f, (t, 0, sp.oo)))       # 1
parcial = sp.simplify(sp.integrate(t * f, (t, 0, R)))
print("media hasta R:", parcial)                    # log(R**2 + 1)/pi
for Rn in (10**3, 10**6):
    print("  R =", Rn, ":", round(float(parcial.subs(R, Rn)), 2))  # 4.4, 8.8
print("media:", sp.integrate(t * f, (t, 0, sp.oo)))            # oo
print("mediana:", sp.tan(sp.pi / 4),
      " p95:", round(float(sp.tan(sp.Rational(19, 40) * sp.pi)), 1))
