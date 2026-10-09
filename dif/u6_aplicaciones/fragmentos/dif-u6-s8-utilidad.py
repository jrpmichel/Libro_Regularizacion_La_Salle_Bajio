# ID: DIF-U6-S8
# Libro: Dif. U6 · subtema 6.8: utilidad máxima, con malla y con la derivada
# Repositorio: dif/u6_aplicaciones/fragmentos/dif-u6-s8-utilidad.py
import numpy as np
import sympy as sp

R = lambda q: 500 * q - 0.5 * q**2          # ingreso, pesos por semana
C = lambda q: 20000 + 140 * q               # costo, pesos por semana
U = lambda q: R(q) - C(q)

q = np.arange(0, 801)                       # malla de unidades enteras
k = np.argmax(U(q))
print("malla: q =", q[k], "U =", U(q[k]))
qs = sp.symbols("q")
Us = sp.nsimplify(U(qs))                    # misma utilidad, exacta
q_opt = sp.solve(sp.diff(Us, qs), qs)[0]    # U' = 360 - q = 0
print("derivada: q =", q_opt, "U =", Us.subs(qs, q_opt),
      "| en q = 500:", Us.subs(qs, 500))
