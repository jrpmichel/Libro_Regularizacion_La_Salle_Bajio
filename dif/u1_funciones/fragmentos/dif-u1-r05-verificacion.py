# ID: DIF-U1-R05
# Libro: Dif. U1 · problema resuelto DIF-U1-05 (verificación)
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-r05-verificacion.py
import numpy as np

V0, R, C = 12.0, 10e3, 47e-6                 # V, ohm, F
tau = R * C
t5 = tau * np.log(20)                        # instante en que v = 0.05 V0
v = lambda t: V0 * np.exp(-t / tau)
print("tau =", round(tau, 4), "s")                                  # 0.47
print("t5 =", round(t5, 4), "s | t5/tau =", round(t5 / tau, 4))     # 1.408 | 2.9957
print(round(v(t5), 4), "V |", round(v(tau), 4), "V")                # 0.6 | 4.4146
print("con log10 en lugar de ln:", round(tau * np.log10(20), 3), "s")   # 0.611, valor equivocado
