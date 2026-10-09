# ID: DIF-U4-R02
# Libro: Dif. U4 · problema resuelto DIF-U4-02 (verificación)
# Repositorio: dif/u4_derivada/fragmentos/dif-u4-r02-verificacion.py
import sympy as sp

g = sp.Rational(981, 100)                       # m/s²
y = sp.symbols("y", positive=True)
h = sp.symbols("h", real=True)
v = sp.sqrt(2 * g * y)
dv = sp.limit((v.subs(y, y + h) - v) / h, h, 0)
print(sp.simplify(dv - g / v))                  # 0: dv/dy = g/v
v2, m = v.subs(y, 2), dv.subs(y, 2)
print(round(float(v2), 4), round(float(m), 4))  # 6.2642 1.566
tan = v2 + m * sp.Rational(2, 10)               # tangente en y = 2.2
exacta = v.subs(y, sp.Rational(22, 10))
print(round(float(tan), 4), round(float(exacta), 4))  # 6.5774 6.5699
