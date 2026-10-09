# ID: DIF-U5-R02
# Libro: Dif. U5 · problema resuelto DIF-U5-02 (verificación)
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-r02-verificacion.py
import sympy as sp

t = sp.symbols("t")                   # meses desde hoy
p = 250 + 4 * t                       # pesos/unidad (lineal, para verificar)
q = 1200 - 15 * t                     # unidades/mes
R = sp.expand(p * q)                  # pesos/mes
print("R(t) =", R)                    # -60*t**2 + 1050*t + 300000
print("p'q =", sp.diff(p, t) * q.subs(t, 0),
      "| pq' =", p.subs(t, 0) * sp.diff(q, t))        # 4800 | -3750
print("R'(0) =", sp.diff(R, t).subs(t, 0))            # 1050
print("R(1) - R(0) =", R.subs(t, 1) - R.subs(t, 0))   # 990 = 1050 - 60
