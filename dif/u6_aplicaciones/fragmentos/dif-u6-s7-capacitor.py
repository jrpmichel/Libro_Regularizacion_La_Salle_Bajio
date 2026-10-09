# ID: DIF-U6-S7
# Libro: Dif. U6 · subtema 6.7: corriente en un capacitor y tensión en un inductor
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-s7-capacitor.py
import sympy as sp

t = sp.symbols("t", real=True)
C, L = 10e-6, 0.02                          # F y H
v = 10 * sp.sin(377 * t)                    # V
i_C = C * sp.diff(v, t)                     # i = C dv/dt
print("i_C(t) =", i_C, "A; pico:", round(10 * 377 * C, 4), "A")
i = 3 * sp.sin(377 * t)                     # A
v_L = L * sp.diff(i, t)                     # v = L di/dt
print("v_L(t) =", v_L, "V; pico:", round(3 * 377 * L, 2), "V")
