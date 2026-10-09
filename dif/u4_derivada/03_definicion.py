# ID: DIF-U4-NB03
# Notebook: dif/u4_derivada.ipynb · sección 4.3 primeros ejemplos
# Repositorio: dif/u4_derivada/03_definicion.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Primeros ejemplos: $x^2$, $x^3$, $1/x$, $\sqrt{x}$
#
# La receta a mano tiene tres pasos: escribir el cociente $\dfrac{f(x+h)-f(x)}{h}$, simplificarlo hasta que $h$ ya no divida (expandir, restar fracciones o multiplicar por el conjugado) y, solo entonces, hacer $h\to0$.
#
# | $f(x)$ | $x^2$ | $x^3$ | $1/x$ | $\sqrt{x}$ |
# |---|---|---|---|---|
# | $f'(x)$ | $2x$ | $3x^2$ | $-1/x^2$ | $\dfrac{1}{2\sqrt{x}}$ ($x>0$) |
#
# La calculadora aplica la definición con `sympy` a la función que escribas y da $f'(x)$ y $f'(a)$. Como comprobación independiente, compara con `sp.diff`, que ya conoce las reglas de la Unidad 5.

# %%
def derivada_por_definicion(f):
    """f'(x) como límite del cociente incremental (variable x)."""
    return sp.simplify(sp.limit((f.subs(x, x + h) - f) / h, h, 0))


def calculadora_definicion(f_txt, a, n_cifras):
    try:
        f = parsear(f_txt)
        d = derivada_por_definicion(f)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print("f(x)  =", f)
    if d.has(sp.nan, sp.zoo, sp.oo, -sp.oo):
        print("f'(x): sympy no obtuvo una fórmula válida para todo x; la función puede no ser derivable")
        print("       en algún punto. Revisa ese punto con las derivadas laterales de la sección 6.")
    else:
        print("f'(x) =", d, "  (por definición)")
        print("comprobación con sp.diff:", "coincide" if sp.simplify(d - sp.diff(f, x)) == 0 else "NO coincide")
    try:
        izq, der, da = derivada_en(f, a)
    except ValueError as err:
        print("En a:", err); return
    print(f"f'({a}) =", (f"{da} ≈ {cifras(da, n_cifras)}  ({n_cifras} cifras significativas)" if da is not None
                         else f"no existe (izquierda: {lado_texto(izq)}; derecha: {lado_texto(der)})"))


widgets.interact(calculadora_definicion,
    f_txt=widgets.Text(value="sqrt(x)", description="f(x) ="),
    a=widgets.FloatText(value=4.0, description="a"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 3 (resultado conocido)
PRUEBAS_3 = [
    ("x² -> 2x", lambda: derivada_por_definicion(x**2) == 2*x),
    ("x³ -> 3x²", lambda: derivada_por_definicion(x**3) == 3*x**2),
    ("1/x -> -1/x²", lambda: sp.simplify(derivada_por_definicion(1/x) + 1/x**2) == 0),
    ("√x -> 1/(2√x)", lambda: sp.simplify(derivada_por_definicion(sp.sqrt(x)) - 1/(2*sp.sqrt(x))) == 0),
    ("3x² - 5x -> 6x - 5", lambda: sp.expand(derivada_por_definicion(3*x**2 - 5*x)) == 6*x - 5),
    ("constante 7 -> 0", lambda: derivada_por_definicion(sp.Integer(7) + 0*x) == 0),
    ("recta 4x - 1 -> 4", lambda: derivada_por_definicion(4*x - 1) == 4),
    ("√x en 4: 1/4", lambda: derivada_en(sp.sqrt(x), 4)[2] == sp.Rational(1, 4)),
]
for nombre, prueba in PRUEBAS_3:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 3).** Calcula a mano, con la definición, la derivada de $f(x)=2x^2+1$. Escribe tu resultado como texto (por ejemplo `"4*x"`) y compáralo con la calculadora.

# %%
mi_derivada = None       # escribe un texto, por ejemplo: "3*x**2"

ref = derivada_por_definicion(2*x**2 + 1)
print("La calculadora da:", ref)
print("falta tu cálculo" if mi_derivada is None else
      ("coincide" if sp.simplify(parsear(mi_derivada) - ref) == 0 else "NO coincide: ¿el 1 se canceló al restar f(x)?"))
