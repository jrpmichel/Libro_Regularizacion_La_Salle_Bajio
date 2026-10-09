# ID: DIF-U6-NB04
# Notebook: dif/u6_aplicaciones.ipynb · sección 6.4 formas indeterminadas y L'Hôpital
# Repositorio: dif/u6_aplicaciones/04_lhopital.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 4. Formas indeterminadas y regla de L'Hôpital
#
# Si $\lim\frac{f}{g}$ es de la forma $\frac00$ o $\frac{\infty}{\infty}$ y $\lim\frac{f'}{g'}$ existe, entonces $\lim\frac fg=\lim\frac{f'}{g'}$. Se derivan **por separado** numerador y denominador, y **antes de cada aplicación** se comprueba que la forma siga siendo indeterminada.
#
# La calculadora aplica la regla paso a paso (hasta cinco veces). Si la forma no es indeterminada, no la aplica y da el límite directo. Para $a=\infty$ escribe `oo` o `∞`. Los resultados son exactos.

# %%
def _forma(num, den, a):
    ln, ld = sp.limit(num, x, a), sp.limit(den, x, a)
    if ln == 0 and ld == 0:
        return "0/0"
    if ln in (sp.oo, -sp.oo) and ld in (sp.oo, -sp.oo):
        return "∞/∞"
    return None


def _punto(a):
    return sp.oo if str(a).strip() in ("oo", "inf", "∞") else sp.nsimplify(a)


def lhopital(num_txt, den_txt, a, max_pasos=5):
    """(límite, pasos) aplicando L'Hôpital mientras la forma sea indeterminada. ValueError si de entrada no lo es."""
    num, den = parsear(num_txt), parsear(den_txt)
    a = _punto(a)
    if _forma(num, den, a) is None:
        raise ValueError(f"la forma no es indeterminada (numerador -> {sp.limit(num, x, a)}, denominador -> "
                         f"{sp.limit(den, x, a)}): no se aplica L'Hôpital; el límite es directo.")
    pasos = []
    while (forma := _forma(num, den, a)) is not None:
        if len(pasos) == max_pasos:
            raise ValueError("después de cinco aplicaciones sigue indeterminado: reescribe el cociente.")
        pasos.append((forma, num, den))
        # se derivan numerador y denominador por separado y se simplifica el nuevo cociente
        num, den = sp.fraction(sp.powsimp(sp.together(sp.diff(num, x) / sp.diff(den, x)), force=True))
    return sp.limit(num / den, x, a), pasos


def calculadora_lhopital(num_txt, den_txt, a):
    try:
        num, den, p = parsear(num_txt), parsear(den_txt), _punto(a)
        if _forma(num, den, p) is None:
            print(f"no es indeterminada (numerador -> {sp.limit(num, x, p)}, denominador -> {sp.limit(den, x, p)}):"
                  f" no se aplica L'Hôpital. Límite directo = {sp.limit(num / den, x, p)}")
            return
        lim, pasos = lhopital(num_txt, den_txt, a)
    except (ValueError, TypeError, sp.SympifyError) as err:
        print("Revisa la entrada:", err); return
    for k, (forma, n_, d_) in enumerate(pasos, 1):
        print(f"paso {k}: ({n_})/({d_}) es {forma} -> derivar por separado")
    print("límite =", lim, "| comprobación con sp.limit:", sp.limit(num / den, x, p))


widgets.interact(calculadora_lhopital,
    num_txt=widgets.Text(value="1 - cos(x)", description="numerador"),
    den_txt=widgets.Text(value="x**2", description="denominador"),
    a=widgets.Text(value="0", description="x ->"));

# %%
# Casos de prueba de la sección 4 (resultado conocido)
PRUEBAS_4 = [
    ("(1 - cos x)/x² en 0: 1/2 con dos pasos", lambda: lhopital("1 - cos(x)", "x**2", 0)[0] == sp.Rational(1, 2)
                                                         and len(lhopital("1 - cos(x)", "x**2", 0)[1]) == 2),
    ("(e^x - 1)/x en 0: 1", lambda: lhopital("exp(x) - 1", "x", 0)[0] == 1),
    ("x³/e^x en ∞: 0 con tres pasos", lambda: lhopital("x**3", "exp(x)", "∞")[0] == 0
                                              and len(lhopital("x**3", "exp(x)", "∞")[1]) == 3),
    ("tan x / x en 0: 1", lambda: lhopital("tan(x)", "x", 0)[0] == 1),
    ("50 ln x / x^0.2 en ∞: 0", lambda: lhopital("50*log(x)", "x**(1/5)", "oo")[0] == 0),
    ("(1 - cos x)/(x² + x) en 0: 0 con un solo paso", lambda: lhopital("1 - cos(x)", "x**2 + x", 0)[0] == 0
                                                         and len(lhopital("1 - cos(x)", "x**2 + x", 0)[1]) == 1),
    ("(x + 1)/(x + 2) en 0 no es indeterminado: se rechaza", lambda: _rechaza(lambda: lhopital("x + 1", "x + 2", 0))),
]
for nombre, prueba in PRUEBAS_4:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 4).** Calcula a mano $\displaystyle\lim_{x\to0}\frac{e^{2x}-1-2x}{x^2}$ y escribe el número. ¿Cuántas veces aplicaste la regla?

# %%
mi_limite = None         # escribe un número

ref, pasos_ref = lhopital("exp(2*x) - 1 - 2*x", "x**2", 0)
if mi_limite is None:                     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    print("La calculadora da:", ref, "con", len(pasos_ref), "aplicaciones")
    print("coincide" if iguales(sp.nsimplify(mi_limite), ref) else "NO coincide: comprueba 0/0 antes de cada aplicación")
