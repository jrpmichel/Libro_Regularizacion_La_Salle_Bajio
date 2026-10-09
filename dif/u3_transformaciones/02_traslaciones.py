# ID: DIF-U3-NB02
# Notebook: dif/u3_transformaciones.ipynb · sección 3.2 traslaciones, reflexiones y paridad
# Repositorio: dif/u3_transformaciones/02_traslaciones.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Traslaciones y reflexiones
#
# - $f(x-h)$ mueve la gráfica $h$ unidades a la **derecha** si $h>0$ (a la izquierda si $h<0$).
# - $f(x)+k$ la mueve $k$ unidades hacia **arriba** si $k>0$.
# - $-f(x)$ la refleja respecto al **eje $x$**; $f(-x)$, respecto al **eje $y$**.
#
# La primera calculadora describe con palabras lo que hace cada cambio y lo comprueba con un punto. La segunda decide si una función es par, impar o ninguna de las dos comparando $f(-x)$ con $f(x)$.

# %%
def describir(h, k, refleja_x, refleja_y):
    """Texto con la secuencia de movimientos de y = s_x·f(s_y·(x - h)) + k."""
    pasos = []
    if refleja_y: pasos.append("reflejar respecto al eje y (cambiar x por -x)")
    if h: pasos.append(f"mover {abs(h)} a la {'derecha' if h > 0 else 'izquierda'}")
    if refleja_x: pasos.append("reflejar respecto al eje x (cambiar y por -y)")
    if k: pasos.append(f"mover {abs(k)} hacia {'arriba' if k > 0 else 'abajo'}")
    return pasos or ["ningún cambio"]


def calculadora_movimientos(base, h, k, refleja_x, refleja_y):
    f, puntos = BASES[base]
    a, b = (-1 if refleja_x else 1), (-1 if refleja_y else 1)
    g = transformar(f, a, b, h, k)
    print("g(x) =", g)
    for i, p in enumerate(describir(h, k, refleja_x, refleja_y), 1):
        print(f"  {i}. {p}")
    p0 = puntos[1]; q = imagen(p0, a, b, h, k)
    print(f"Comprobación con un punto: {p0} de f  ->  {q} de g;  g({q[0]}) = {sp.simplify(g.subs(x, q[0]))}")


widgets.interact(calculadora_movimientos,
    base=widgets.Dropdown(options=list(BASES), value="x²", description="base"),
    h=widgets.IntSlider(value=3, min=-5, max=5, description="h"),
    k=widgets.IntSlider(value=-1, min=-5, max=5, description="k"),
    refleja_x=widgets.Checkbox(value=False, description="reflejar en eje x"),
    refleja_y=widgets.Checkbox(value=False, description="reflejar en eje y"));

# %%
def paridad(texto):
    """'par', 'impar' o 'ninguna', comparando f(-x) con f(x) de forma simbólica."""
    f = parsear(texto)
    fm = f.subs(x, -x)
    if sp.simplify(fm - f) == 0:
        return "par"
    if sp.simplify(fm + f) == 0:
        return "impar"
    return "ninguna"


def calculadora_paridad(texto):
    try:
        tipo = paridad(texto)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"f(x) = {texto}:  {tipo}")
    print("par: simétrica respecto al eje y;  impar: simétrica respecto al origen.")


widgets.interact(calculadora_paridad, texto=widgets.Text(value="x*cos(x)", description="f(x) ="));

# %%
# Casos de prueba de la sección 2 (resultado conocido)
PRUEBAS_2 = [
    ("x³ es impar, cos x es par, x² + x no es ninguna",
     lambda: (paridad("x**3"), paridad("cos(x)"), paridad("x**2 + x")) == ("impar", "par", "ninguna")),
    ("x cos x es impar y e^(-x²) es par", lambda: (paridad("x*cos(x)"), paridad("exp(-x**2)")) == ("impar", "par")),
    ("(x+2)² - 1 tiene vértice (-2, -1)", lambda: imagen((0, 0), 1, 1, -2, -1) == (-2, -1)),
    ("√(-x) tiene dominio (-∞, 0]", lambda: continuous_domain(sp.sqrt(-x), x, sp.S.Reals) == sp.Interval(-sp.oo, 0)),
    ("e^(-x) en x = 1 vale lo que e^x en x = -1", lambda: transformar(sp.exp(x), 1, -1, 0, 0).subs(x, 1) == sp.exp(-1)),
    ("h = 3 se describe como 'mover 3 a la derecha'", lambda: describir(3, 0, False, False) == ["mover 3 a la derecha"]),
    ("una expresión con otra letra se rechaza", lambda: _rechaza(lambda: paridad("t**2"))),
]
for nombre, prueba in PRUEBAS_2:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 2).** Escribe a mano el vértice de $y=(x+3)^2+1$ y di hacia dónde se movió la parábola $y=x^2$. Compara con la calculadora (base $x^2$, $h=-3$, $k=1$).

# %%
mi_vertice = None       # escribe una pareja, por ejemplo: (0, 0)

ref = imagen((0, 0), 1, 1, -3, 1)
print("La calculadora da:", ref)
print("falta tu cálculo" if mi_vertice is None else
      ("coincide" if tuple(mi_vertice) == ref else "NO coincide: ¿qué valor de x hace cero a x + 3?"))
