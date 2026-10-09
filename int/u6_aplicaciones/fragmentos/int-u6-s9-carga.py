# ID: INT-U6-S9
# Libro: Int. U6 · subtema 6.9: carga total y valor eficaz de un PWM
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-s9-carga.py
import sympy as sp

t = sp.symbols("t", nonnegative=True)
i = sp.Rational(1, 2) * sp.exp(-50 * t)        # A, descarga de INT-U4-02
print("carga total:", sp.integrate(i, (t, 0, sp.oo)), "C")      # 1/100
V, D, T = sp.symbols("V D T", positive=True)
media = sp.integrate(V, (t, 0, D * T)) / T     # PWM: vale V durante D*T
eficaz = sp.sqrt(sp.integrate(V**2, (t, 0, D * T)) / T)
print("PWM: media", media, "; eficaz", eficaz)  # D*V; sqrt(D)*V
