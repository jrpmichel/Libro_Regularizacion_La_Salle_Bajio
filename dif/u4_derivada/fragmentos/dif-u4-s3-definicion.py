# ID: DIF-U4-S3
# Libro: Dif. U4 · subtema 4.3: derivada por definición con sympy
# Repositorio: dif/u4_derivada/fragmentos/dif-u4-s3-definicion.py
import sympy as sp

x = sp.symbols("x", positive=True)
h = sp.symbols("h", real=True)
for f in [x**2, x**3, 1/x]:
    cociente = sp.cancel((f.subs(x, x + h) - f) / h)   # h ya no divide
    print(f"{str(f):4}  cociente = {cociente}   h -> 0: "
          f"{cociente.subs(h, 0)}")

f = sp.sqrt(x)                                    # con raíz: conjugado
conj = sp.sqrt(x + h) + sp.sqrt(x)
num = sp.expand((f.subs(x, x + h) - f) * conj)    # (x + h) - x = h
cociente = num / (h * conj)
print("sqrt(x)  cociente =", cociente, "  h -> 0:", cociente.subs(h, 0))
