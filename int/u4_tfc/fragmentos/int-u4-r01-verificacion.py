# ID: INT-U4-R01
# Libro: Int. U4 · problema resuelto INT-U4-01 (verificación)
# Repositorio: int/u4_tfc/fragmentos/int-u4-r01-verificacion.py
import sympy as sp

t = sp.symbols("t", real=True)
Q = sp.Piecewise((3 * t, t <= 3), (18 - 3 * t, True))   # miles de m^3/h
V0 = 120                                                 # miles de m^3
V = lambda T: V0 + sp.integrate(Q, (t, 0, T))
print("Q = 0 en t =", sp.solve(18 - 3 * t, t))          # [6]
for T in (3, 6, 8):
    print(f"V({T}) =", float(V(T)))                      # 133.5, 147, 141
print("área de 6 a 8 h:", sp.integrate(Q, (t, 6, 8)))   # -6
