# ID: DIF-U1-NB03
# Notebook: dif/u1_funciones.ipynb · sección 1.3 dominio y rango
# Repositorio: dif/u1_funciones/03_dominio_rango.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Dominio y rango
#
# La calculadora devuelve el dominio y el rango con `sympy` y además barre $[-10,10]$ buscando puntos donde la regla no produce un número real. El barrido es una comprobación aproximada: confirma al dominio, no lo demuestra.
#
# Condiciones que aplica la teoría: denominador distinto de cero y radicando de raíz de índice par mayor o igual que cero, sobre la expresión **sin simplificar**.

# %%
def dominio_rango(texto):
    f = parsear(texto)
    try:
        d = continuous_domain(f, x, sp.S.Reals)
    except NotImplementedError:
        d = None
    try:
        r = function_range(f, x, sp.S.Reals)
    except NotImplementedError:
        r = None
    return f, d, r


def calculadora_dominio(texto):
    try:
        f, d, r = dominio_rango(texto)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"f(x) = {texto}")
    print("Dominio:", "sympy no pudo determinarlo" if d is None else conjunto_a_texto(d))
    print("Rango  :", "sympy no pudo determinarlo" if r is None else conjunto_a_texto(r))
    fn = sp.lambdify(x, f, "numpy")
    xs = np.linspace(-10, 10, 2001)
    with np.errstate(all="ignore"):
        ys = fn(xs)
    ys = np.full_like(xs, float(ys)) if np.ndim(ys) == 0 else np.asarray(ys, dtype=float)
    malos = xs[~np.isfinite(ys)]
    print("Barrido en [-10, 10]:",
          "f produce valores reales en todos los puntos revisados." if malos.size == 0
          else f"f falla en {malos.size} de {xs.size} puntos revisados (de x ≈ {cifras(malos.min(), 4)} a x ≈ {cifras(malos.max(), 4)}).")
    plt.figure(figsize=(4.5, 3)); plt.plot(xs, np.where(np.isfinite(ys) & (np.abs(ys) < 20), ys, np.nan), "b-")
    plt.grid(True); plt.xlabel("x"); plt.ylabel("y"); plt.ylim(-6, 6); plt.show()


widgets.interact(calculadora_dominio, texto=widgets.Text(
    value="sqrt(x + 3)/(x**2 - 4)", description="f(x) =", layout=widgets.Layout(width="420px")));

# %%
def en_texto(f):
    _, d, r = dominio_rango(f)
    return conjunto_a_texto(d), conjunto_a_texto(r)

PRUEBAS_3 = [
    ("√(x+3): dominio [-3, ∞), rango [0, ∞)", lambda: en_texto("sqrt(x+3)") == ("[-3, ∞)", "[0, ∞)")),
    ("1/(x-2): dominio y rango sin 2 y sin 0",  lambda: en_texto("1/(x-2)") == ("(-∞, 2) ∪ (2, ∞)", "(-∞, 0) ∪ (0, ∞)")),
    ("√(9-x²): dominio [-3,3], rango [0,3]",    lambda: en_texto("sqrt(9-x**2)") == ("[-3, 3]", "[0, 3]")),
    ("√(x+3)/(x²-4): dominio [-3,-2) ∪ (-2,2) ∪ (2,∞)", lambda: en_texto("sqrt(x+3)/(x**2-4)")[0] == "[-3, -2) ∪ (-2, 2) ∪ (2, ∞)"),
    ("(x²-4)/(x-2) conserva el hueco en x=2",    lambda: en_texto("(x**2-4)/(x-2)")[0] == "(-∞, 2) ∪ (2, ∞)"),
]
for nombre, prueba in PRUEBAS_3:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 3).** Determina a mano el dominio de $f(x)=\dfrac{\sqrt{x-1}}{x-4}$. Escríbelo como intervalos, ejecuta la calculadora con esa expresión y compara. Anota: ¿incluiste el extremo $1$? ¿excluiste el $4$?

# %%
mi_dominio = ""        # por ejemplo: "[1, 4) ∪ (4, ∞)"
_, d, _ = dominio_rango("sqrt(x - 1)/(x - 4)")
print("La calculadora da:", conjunto_a_texto(d))
print("Tu resultado     :", mi_dominio or "(falta tu cálculo a mano)")
