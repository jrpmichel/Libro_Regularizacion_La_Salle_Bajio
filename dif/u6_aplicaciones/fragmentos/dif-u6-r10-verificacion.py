# ID: DIF-U6-R10
# Libro: Dif. U6 · problema resuelto DIF-U6-10 (verificación)
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-r10-verificacion.py
import sympy as sp

Q, D, S, H = sp.symbols("Q D S H", positive=True)
CT = D * S / Q + H * Q / 2                        # pesos por año
Qs = sp.solve(sp.diff(CT, Q), Q)[0]
print("Q* =", Qs, "| CT'' =", sp.diff(CT, Q, 2))
base = {D: 12000, H: 30}
A, B = {**base, S: 800}, {**base, S: 1100}        # proveedores A y B
print("A: Q* =", Qs.subs(A), "| CT(600) =", CT.subs(A).subs(Q, 600))
qB = Qs.subs(B)
print("B: Q* =", round(float(qB), 1),
      "| CT =", round(float(CT.subs(B).subs(Q, qB)), 1))
S_lim = sp.solve(CT.subs(base).subs(Q, Qs.subs(base)) - 25000, S)[0]
print("B conviene si S <", round(float(S_lim), 1))
# A: 800 > 600 y 25000 | B: 938.1 y 28142.5 | S < 868.1
