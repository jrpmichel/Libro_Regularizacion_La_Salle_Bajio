# ID: DIF-U1-NB04
# Notebook: dif/u1_funciones.ipynb · sección 1.4 álgebra de apoyo
# Repositorio: dif/u1_funciones/04_algebra.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 4. Álgebra de apoyo sobre la marcha
#
# La calculadora expande $f(x+h)$, simplifica el cociente de diferencias $\dfrac{f(x+h)-f(x)}{h}$ y comprueba que $f(x+h)$ **no** es $f(x)+h$. Después muestra por qué conviene simplificar antes de evaluar: con $h$ muy pequeño, restar dos números casi iguales destruye cifras significativas.

# %%
def cociente_diferencias(texto):
    f = parsear(texto)
    fxh = sp.expand(f.subs(x, x + h))
    cociente = sp.simplify((fxh - f) / h)
    es_mas_h = sp.simplify(fxh - (f + h)) == 0
    return f, fxh, cociente, es_mas_h


def calculadora_cociente(texto, x0, exponente_h):
    try:
        f, fxh, cociente, es_mas_h = cociente_diferencias(texto)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"f(x)               = {f}")
    print(f"f(x+h) expandida   = {fxh}")
    print(f"[f(x+h)-f(x)]/h    = {cociente}      (válido para h ≠ 0)")
    print(f"¿f(x+h) = f(x)+h?  = {es_mas_h}")
    fn = sp.lambdify(x, f, "numpy"); simp = sp.lambdify((x, h), cociente, "numpy")
    print(f"\nEn x = {x0}: cociente directo contra cociente simplificado")
    print(f"{'h':>10} | {'directo':>14} | {'simplificado':>14}")
    for e in sorted({exponente_h, 1, 4, 8, 12, 16}):
        hh = 10.0**(-e)
        directo = (fn(x0 + hh) - fn(x0)) / hh
        print(f"{hh:>10.0e} | {directo:>14.8g} | {float(simp(x0, hh)):>14.8g}")


widgets.interact(calculadora_cociente,
    texto=widgets.Text(value="x**2 - 3*x", description="f(x) ="),
    x0=widgets.FloatText(value=4.0, description="x"),
    exponente_h=widgets.IntSlider(value=6, min=1, max=16, description="h = 10^-k"));

# %%
# Precaución de herramienta: la jerarquía de operaciones también se aplica en Python
for v in (4, 9):
    print(f"x = {v}:  x**1/2 = {v**1/2:<5}  x**(1/2) = {v**(1/2)}")
# En x = 4 los dos coinciden por casualidad (4/2 = 2 = raíz de 4); en x = 9 ya no.

# %%
PRUEBAS_4 = [
    ("x² -> 2x + h",            lambda: sp.simplify(cociente_diferencias("x**2")[2] - (2*x + h)) == 0),
    ("x³ -> 3x² + 3xh + h²",    lambda: sp.simplify(cociente_diferencias("x**3")[2] - (3*x**2 + 3*x*h + h**2)) == 0),
    ("2x + 1 -> 2",             lambda: sp.simplify(cociente_diferencias("2*x + 1")[2] - 2) == 0),
    ("x²-3x -> 2x + h - 3",     lambda: sp.simplify(cociente_diferencias("x**2 - 3*x")[2] - (2*x + h - 3)) == 0),
    ("f(x+h) no es f(x)+h para x²-3x", lambda: cociente_diferencias("x**2 - 3*x")[3] is False),
]
for nombre, prueba in PRUEBAS_4:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 4).** Calcula a mano $\dfrac{f(x+h)-f(x)}{h}$ para $f(x)=x^2-3x$. Compara con la calculadora. Anota: ¿cancelaste $h$ como factor? ¿qué condición sobre $h$ necesitas para poder hacerlo?

# %%
mi_resultado = ""      # por ejemplo: "2*x + h - 3"
comparar(mi_resultado, cociente_diferencias("x**2 - 3*x")[2])
