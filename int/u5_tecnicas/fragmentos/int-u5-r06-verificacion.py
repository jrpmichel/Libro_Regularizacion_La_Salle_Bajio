# ID: INT-U5-R06
# Libro: Int. U5 · problema resuelto INT-U5-06 (verificación)
# Repositorio: int/u5_tecnicas/fragmentos/int-u5-r06-verificacion.py
import sympy as sp

t = sp.symbols("t")
Vp = sp.symbols("V_p", positive=True)
T = sp.Rational(1, 60)                               # 60 Hz
v = Vp * sp.sin(120 * sp.pi * t)
Vrms = sp.sqrt(sp.integrate(v**2, (t, 0, T)) / T)
print("V_rms =", sp.simplify(Vrms))                  # sqrt(2)*V_p/2
print("V_p para 127 V:", round(127 * 2**0.5, 1), "V")    # 179.6
print("P =", round(127**2 / 20, 2), "W")                 # 806.45
