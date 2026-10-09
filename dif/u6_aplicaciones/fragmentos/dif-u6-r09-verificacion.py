# ID: DIF-U6-R09
# Libro: Dif. U6 · problema resuelto DIF-U6-09 (verificación)
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-r09-verificacion.py
import sympy as sp

q = sp.symbols("q", positive=True)        # unidades por semana
R = q * (500 - q / 2)                     # ingreso, pesos por semana
U = R - (20000 + 140 * q)                 # utilidad
qR = sp.solve(sp.diff(R, q), q)[0]
qU = sp.solve(sp.diff(U, q), q)[0]
print("máx. ingreso: q =", qR, " U =", U.subs(q, qR))       # 500, 35000
print("máx. utilidad: q =", qU, " U =", U.subs(q, qU))      # 360, 44800
print("precio:", 500 - qU / 2, "| U'' =", sp.diff(U, q, 2))  # 320, -1
