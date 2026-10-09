# ID: INT-U3-L3
# Libro: Int. U3 · Laboratorio del Error 3.3 (prueba con código)
# Repositorio: int/u3_antiderivada/fragmentos/int-u3-l3-prueba.py
import sympy as sp

t, C1, C2 = sp.symbols("t C1 C2")
a = sp.Rational(6, 5) - sp.Rational(6, 25) * t       # 1.2 - 0.24 t, m/s^2
v = sp.integrate(a, t) + C1
s = sp.integrate(v, t) + C2
s_ia = s.subs({C1: 0, C2: 0})                         # la IA fija C1 = 0
v_ok = v.subs(C1, sp.Rational(1, 2))                  # v(0) = 0.5 m/s
s_ok = s.subs({C1: sp.Rational(1, 2), C2: 0})         # s(0) = 0
print("IA:   s(5) =", s_ia.subs(t, 5), "m")
print("bien: s(5) =", s_ok.subs(t, 5), "m   v(5) =", v_ok.subs(t, 5), "m/s")
print("revisión: v(0) =", v_ok.subs(t, 0), "  s(0) =", s_ok.subs(t, 0))
