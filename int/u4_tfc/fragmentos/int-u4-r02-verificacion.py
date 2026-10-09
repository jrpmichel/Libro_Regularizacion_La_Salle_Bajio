# ID: INT-U4-R02
# Libro: Int. U4 · problema resuelto INT-U4-02 (verificación)
# Repositorio: int/u4_tfc/fragmentos/int-u4-r02-verificacion.py
import sympy as sp

t = sp.symbols("t", nonnegative=True)
R, Cc, V0 = 100, sp.Rational(200, 10**6), 50     # ohm, F, V
tau = R * Cc                                      # 0.02 s
i = V0 / R * sp.exp(-t / tau)                     # 0.5 e^(-50 t) A
E = sp.integrate(R * i**2, (t, 0, sp.Rational(5, 100)))
print("E =", sp.simplify(E), "=", round(float(E), 4), "J")   # 0.2483
print("guardada:", Cc * V0**2 / 2, "J")                      # 1/4
print("fracción:", round(float(E / (Cc * V0**2 / 2)), 4))    # 0.9933
