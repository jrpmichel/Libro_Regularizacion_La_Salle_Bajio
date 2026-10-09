# ID: INT-U6-R02
# Libro: Int. U6 · problema resuelto INT-U6-02 (verificación)
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-r02-verificacion.py
import sympy as sp

x = sp.symbols("x", real=True)
r = 5 - x**2 / 2000                        # radio exterior (cm)
ri = sp.Rational(42, 10)                    # radio interior (cm)
V = sp.pi * sp.integrate(r**2 - ri**2, (x, -20, 20))
print("V =", V, "=", round(float(V), 1), "cm^3")         # 842.1
print("masa:", round(2.70 * float(V)), "g")              # 2274
for R in (5, sp.Rational(48, 10)):                       # cotas
    print("cilindro de radio", R, ":",
          round(float(sp.pi * (R**2 - ri**2) * 40), 1), "cm^3")
