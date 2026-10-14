# ID: ALG-U6-R05
# Libro: Alg. U6 · problema resuelto ALG-U6-05 (verificacion)
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-r05-verificacion.py
import sympy as sp

h = sp.symbols("h", positive=True)
L = sp.sqrt(9 + h**2)
AB, BC, AD, CD, BD, Ax, Ay, Cy = X = sp.symbols("AB BC AD CD BD Ax Ay Cy")
c, s = 3 / L, h / L
ecs = [AB + c * AD + Ax, s * AD + Ay, -AB + BC, BD - 12,
       -BC - c * CD, s * CD + Cy, -c * AD + c * CD + 3, -s * AD - s * CD - BD]
f = {k: sp.simplify(v) for k, v in sp.solve(ecs, X, dict=True)[0].items()}
print("AB =", f[AB], "; AD =", f[AD], "; CD =", f[CD])
# para 0 < h < 12: AB, BC, BD en tension; AD y CD en compresion
S = sp.simplify(3 * f[AB] + 3 * f[BC] - L * f[AD] - L * f[CD] + h * f[BD])
print("S(h) =", sp.expand(S))
print("dS/dh = 0 en h =", sp.solve(sp.diff(S, h), h),
      "; S(3) =", S.subs(h, 3), "; S(4) =", S.subs(h, 4))
comp = -f[CD]                                   # la barra mas comprimida
for hh in (3, 4):
    print(f"h = {hh}: CD = {float(comp.subs(h, hh)):.3f} kN,",
          f"P_adm = {float(300 / L.subs(h, hh)**2):.2f} kN")
lim = [sp.nsolve(comp * L**2 - 300, h, g) for g in (0.6, 4.6)]
print("pandeo bien para h entre", [round(float(v), 3) for v in lim], "m")
print(f"acero: {float(S.subs(h, 3)) * 1e3 / 150e6 * 7850:.1f} kg")
