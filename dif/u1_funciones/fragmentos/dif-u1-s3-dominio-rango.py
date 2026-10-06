# ID: DIF-U1-S3
# Libro: Dif. U1 · subtema 1.3: dominio y rango con sympy
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-s3-dominio-rango.py
import sympy as sp
from sympy.calculus.util import continuous_domain, function_range

x = sp.symbols("x", real=True)
for f in (sp.sqrt(x + 3), 1/(x - 2)):
    print(f, "| dominio:", continuous_domain(f, x, sp.S.Reals),
          "| rango:", function_range(f, x, sp.S.Reals))
