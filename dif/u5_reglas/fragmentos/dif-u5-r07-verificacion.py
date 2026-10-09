# ID: DIF-U5-R07
# Libro: Dif. U5 · problema resuelto DIF-U5-07 (verificación)
# Repositorio: dif/u5_reglas/fragmentos/dif-u5-r07-verificacion.py
import sympy as sp

d, L, rho, ed, eL = sp.symbols("d L rho e_d e_L", positive=True)
m = rho * sp.pi * d**2 / 4 * L
rel = sp.diff(sp.log(m), d) * ed + sp.diff(sp.log(m), L) * eL
print("dm/m =", sp.simplify(rel))                 # 2*e_d/d + e_L/L
casos = {"actual": (0.02, 0.5), "micrometro": (0.002, 0.5),
         "calibrador": (0.02, 0.02)}                # mm
for nombre, (a, b) in casos.items():
    v = rel.subs({d: 25, L: 120, ed: a, eL: b})
    print(f"{nombre:10}: {100 * float(v):.3f} %")  # 0.577, 0.433, 0.177
cambio = m.subs({d: 25.02, L: 120.5}) / m.subs({d: 25, L: 120}) - 1
print("directo   :", round(100 * float(cambio), 3), "%")   # 0.577
