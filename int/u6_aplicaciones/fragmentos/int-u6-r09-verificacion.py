# ID: INT-U6-R09
# Libro: Int. U6 · problema resuelto INT-U6-09 (verificación)
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-r09-verificacion.py
import sympy as sp

t = sp.symbols("t", nonnegative=True)
q_cal = sp.integrate(3 + 9 * sp.exp(-t / sp.Rational(2, 5)),
                     (t, 0, sp.Rational(3, 2)))            # mC
q_act = float(q_cal) + 40 * 0.12                           # + transmisión
print(f"calentamiento {float(q_cal):.3f} mC; activo {q_act:.3f} mC")
Q = 2600 * 0.8 * 3.6                                       # C útiles
I_max = Q / (2 * 365 * 24 * 3600)                          # A, 2 años
T = (q_act * 1e-3 - 6e-6 * 1.62) / (I_max - 6e-6)
print(f"corriente máxima {I_max * 1e6:.1f} uA; periodo mínimo {T:.1f} s")
for P in (120, 300):
    I = (q_act * 1e-3 + 6e-6 * (P - 1.62)) / P
    print(f"periodo {P} s: {I * 1e6:.1f} uA, {Q / I / 3.1536e7:.2f} años")
