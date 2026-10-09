# ID: INT-U3-R03
# Libro: Int. U3 · problema resuelto INT-U3-03 (verificación)
# Repositorio: int/u3_antiderivada/fragmentos/int-u3-r03-verificacion.py
import sympy as sp

x = sp.symbols("x", positive=True)            # (a) vale para x > 0
t = sp.symbols("t", real=True)                # (b) vale a cada lado del 0
fa = 6*sp.sqrt(x) - 4/x**2 + 2*sp.exp(-x/2)
Fa = 4*x**sp.Rational(3, 2) + 4/x - 4*sp.exp(-x/2)
fb = (t**2 - 1)/t + 6*sp.cos(3*t)
Fb = t**2/2 - sp.log(sp.Abs(t)) + 2*sp.sin(3*t)
print("(a) F' - f =", sp.simplify(sp.diff(Fa, x) - fa))          # 0
for c in (-1, 2):                         # (b) en t < 0 y en t > 0
    d = sp.diff(Fb, t).subs(t, c) - fb.subs(t, c)
    print(f"(b) F' - f en t = {c}:", sp.simplify(d))              # 0
print("integrando (a) en x = 1:", round(float(fa.subs(x, 1)), 3))  # 3.213
