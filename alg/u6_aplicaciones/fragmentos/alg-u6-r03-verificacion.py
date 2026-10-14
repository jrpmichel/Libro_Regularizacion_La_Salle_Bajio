# ID: ALG-U6-R03
# Libro: Alg. U6 · problema resuelto ALG-U6-03 (verificacion)
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-r03-verificacion.py
import sympy as sp

pct = lambda *v: [sp.Rational(x, 100) for x in v]
Cu, Zn, Pb = pct(80, 60, 65), pct(18, 35, 25), pct(2, 5, 10)
S = sp.Matrix([[1, 1, 1], Cu, Zn])          # columnas: chatarras 1, 2, 3
for cu, zn in ((70, 25), (70, 28), (75, 22)):
    x = S.solve(sp.Matrix([1000, 10 * cu, 10 * zn]))   # kg para 1000 kg
    plomo = sum(p * xi for p, xi in zip(Pb, x))
    print(f"{cu}/{zn}: x =", [round(float(v), 1) for v in x],
          "factible" if min(x) >= 0 else "NO factible",
          f"| plomo {float(plomo) / 10:.1f} %")
print("columnas suman 1:", [a + b + c for a, b, c in zip(Cu, Zn, Pb)])
