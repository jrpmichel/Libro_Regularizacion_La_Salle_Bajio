# ID: DIF-U4-R03
# Libro: Dif. U4 · problema resuelto DIF-U4-03 (verificación)
# Repositorio: dif/u4_derivada/fragmentos/dif-u4-r03-verificacion.py
import sympy as sp

t, h = sp.symbols("t h", real=True)
s = sp.Rational(6, 10) * t**2 - sp.Rational(1, 10) * t**3   # m, s
cociente = sp.expand((s.subs(t, 2 + h) - s.subs(t, 2)) / h)
print(cociente)                             # 6/5 - h**2/10
print(sp.limit(cociente, h, 0))             # 6/5 = 1.2 m/s
print(cociente.subs(h, 1))                  # 11/10: dato del compañero
print(sp.solve(sp.Eq(h**2 / 10, sp.Rational(1, 100)), h))  # +-sqrt(10)/10
print(cociente.subs(h, sp.Rational(1, 4)))  # 191/160 = 1.19375: cada 0.25 s
