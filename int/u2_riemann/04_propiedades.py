# ID: INT-U2-NB04
# Notebook: int/u2_riemann.ipynb · sección 2.4 propiedades de la integral definida
# Repositorio: int/u2_riemann/04_propiedades.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 4. Propiedades básicas de la integral definida
#
# * Convenios: $\int_a^a f=0$ y $\int_b^a f=-\int_a^b f$.
# * Linealidad: $\int(\alpha f+\beta g)=\alpha\int f+\beta\int g$. Aditividad: $\int_a^b=\int_a^c+\int_c^b$ (aunque $c$ quede fuera).
# * Cotas: si $m\le f\le M$ en $[a,b]$, entonces $m(b-a)\le\int_a^b f\le M(b-a)$.
# * Simetría en $[-a,a]$: par, $2\int_0^a f$; impar, $0$.
#
# La calculadora comprueba cada propiedad con `sympy` (exacto) y con una suma de punto medio (numérico). Que coincidan las dos vías es la comprobación.

# %%
def integral_exacta(f, a, b):
    return sp.integrate(f, (x, sp.nsimplify(a), sp.nsimplify(b)))


def integral_numerica(f, a, b, m=20000):
    """Suma de punto medio con m franjas; acepta a > b (Δx negativo)."""
    F = numerica(f)
    xs = np.linspace(a, b, m + 1)
    return float(F((xs[:-1] + xs[1:]) / 2).sum() * (b - a) / m)


def extremos(f, a, b):
    """(m, M) de f en [a, b] con los extremos del intervalo y los puntos donde f' = 0."""
    cand = [sp.nsimplify(a), sp.nsimplify(b)]
    crit = sp.solveset(sp.diff(f, x), x, domain=sp.Interval.open(sp.nsimplify(a), sp.nsimplify(b)))
    if isinstance(crit, sp.FiniteSet):
        cand += list(crit)
    elif crit is not sp.S.EmptySet:
        cand += list(np.linspace(float(a), float(b), 401))
    vals = [valor_real(f, c) for c in cand]
    return min(vals), max(vals)


def revisar_propiedades(f, g, a, b, c, alfa, beta):
    """Diccionario {propiedad: (lado izquierdo, lado derecho)} evaluado con sympy."""
    I = lambda h, p, q: integral_exacta(h, p, q)
    r = {"linealidad": (I(alfa * f + beta * g, a, b), alfa * I(f, a, b) + beta * I(g, a, b)),
         "aditividad": (I(f, a, b), I(f, a, c) + I(f, c, b)),
         "límites al revés": (I(f, b, a), -I(f, a, b))}
    if a < b:
        m_, M_ = extremos(f, a, b)
        r["cotas"] = (m_ * (b - a), float(I(f, a, b)), M_ * (b - a))
    return r


def calculadora_propiedades(f_txt, g_txt, a, b, c, alfa, beta, n_cifras):
    try:
        f, g = parsear(f_txt), parsear(g_txt)
        r = revisar_propiedades(f, g, a, b, c, alfa, beta)
    except (ValueError, TypeError) as err:
        print("Revisa la entrada:", err); return
    for nombre, valores in r.items():
        if nombre == "cotas":
            print(f"  cotas: {cifras(valores[0], n_cifras)} <= {cifras(valores[1], n_cifras)} <= {cifras(valores[2], n_cifras)}")
        else:
            izq, der = valores
            ok = "se cumple" if sp.simplify(izq - der) == 0 else "NO se cumple"
            print(f"  {nombre}: {cifras(izq, n_cifras)} y {cifras(der, n_cifras)} -> {ok}")
    print(f"  comprobación numérica de ∫ f de a a b (punto medio): {cifras(integral_numerica(f, a, b), n_cifras)}")
    print(f"  ({n_cifras} cifras significativas, redondeado)")


widgets.interact(calculadora_propiedades,
    f_txt=widgets.Text(value="x**2", description="f(x) ="), g_txt=widgets.Text(value="x", description="g(x) ="),
    a=widgets.FloatText(value=0, description="a"), b=widgets.FloatText(value=2, description="b"),
    c=widgets.FloatText(value=3, description="c"),
    alfa=widgets.FloatText(value=3, description="α"), beta=widgets.FloatText(value=-2, description="β"),
    n_cifras=widgets.IntSlider(value=5, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 4 (resultado conocido)
PRUEBAS_4 = [
    ("linealidad: ∫(3x^2 - 2x + 5) en [0, 1] = 5", lambda: integral_exacta(3 * x**2 - 2 * x + 5, 0, 1) == 5),
    ("aditividad con c fuera: ∫_1^2 x^2 = ∫_1^3 + ∫_3^2 = 7/3",
     lambda: integral_exacta(x**2, 1, 3) + integral_exacta(x**2, 3, 2) == sp.Rational(7, 3)),
    ("límites al revés: ∫_2^1 x^2 = -7/3, también con sumas",
     lambda: integral_exacta(x**2, 2, 1) == -sp.Rational(7, 3) and cerca(integral_numerica(x**2, 2, 1), -7 / 3, 1e-6)),
    ("cotas de e^(-x^2) en [0, 1]: e^-1 y 1", lambda: cerca(extremos(sp.exp(-x**2), 0, 1)[0], math.exp(-1))
                                                   and cerca(extremos(sp.exp(-x**2), 0, 1)[1], 1)),
    ("simetría: ∫x^3 en [-2, 2] = 0 y ∫x^2 en [-2, 2] = 16/3",
     lambda: integral_exacta(x**3, -2, 2) == 0 and integral_exacta(x**2, -2, 2) == sp.Rational(16, 3)),
    ("la cota de x^2 + 1 en [0, 1] descarta 1/3", lambda: revisar_propiedades(x**2 + 1, x, 0, 1, 2, 1, 1)["cotas"][0] > 1 / 3),
]
for nombre, prueba in PRUEBAS_4:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 4).** Calcula a mano $\int_{-1}^{1}\big(x^3+2\big)dx$ con la simetría y la integral de una constante.

# %%
mi_integral = None     # escribe un número

if mi_integral is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    ref = integral_exacta(x**3 + 2, -1, 1)
    print("La calculadora da:", ref)
    print("coinciden" if cerca(mi_integral, ref, 1e-9) else "NO coinciden: x^3 es impar y la constante 2 aporta 2·(1 - (-1))")
