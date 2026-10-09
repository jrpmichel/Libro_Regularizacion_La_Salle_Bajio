# ID: DIF-U3-S2
# Libro: Dif. U3 · subtema 3.2: traslaciones, reflexiones y paridad
# Repositorio: dif/u3_transformaciones/fragmentos/dif-u3-s2-traslaciones-reflexiones.py
import math

f = lambda x: x**2
g = lambda x: f(x - 3)                    # ¿derecha o izquierda?
for x in (0, 3, 6):
    print(f"g({x}) = f({x - 3}) = {g(x)}")    # f(0) = 0 aparece en x = 3

print([math.sqrt(-x) for x in (-4, -1, 0)])  # sqrt(-x): si x <= 0
print([-math.sqrt(x) for x in (1, 4, 9)])    # -sqrt(x): si x >= 0

# paridad: f(-x) = f(x) es par; f(-x) = -f(x) es impar
pruebas = (0.5, 1.0, 2.0)
for nombre, p in [("x^3", lambda x: x**3), ("cos x", math.cos),
                  ("x^2 + x", lambda x: x**2 + x)]:
    par = all(math.isclose(p(-x), p(x)) for x in pruebas)
    impar = all(math.isclose(p(-x), -p(x)) for x in pruebas)
    print(nombre, "->", "par" if par else "impar" if impar else "ninguna")
