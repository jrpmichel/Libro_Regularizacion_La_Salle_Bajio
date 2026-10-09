# ID: INT-U4-NB02
# Notebook: int/u4_tfc.ipynb · sección 4.2 segunda parte: evaluación con la antiderivada
# Repositorio: int/u4_tfc/02_segunda.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Segunda parte: evaluación con la antiderivada
#
# $$\int_a^b f(x)\,dx=F(b)-F(a)=\Big[F(x)\Big]_a^b,\qquad\text{si } f \text{ es continua en } [a,b] \text{ y } F'=f.$$
#
# La calculadora revisa primero la continuidad de $f$ en $[a,b]$; si no se cumple, **no** aplica el teorema y explica por qué. Después busca una antiderivada, muestra la resta y la compara con una suma de punto medio con 1000 franjas.

# %%
def evaluar(f, a, b):
    """(F, F(b) - F(a)) con la segunda parte del teorema; ValueError si no aplica."""
    if a == b:
        return sp.Integer(0) * x, sp.Integer(0)
    if not continua_en(f, a, b):
        raise ValueError("f no es continua en todo el intervalo [a, b]: el teorema fundamental no aplica "
                         "y F(b) - F(a) puede dar un número sin sentido.")
    an, bn = sp.nsimplify(a), sp.nsimplify(b)
    if f.has(sp.Abs):                          # |u|: se separa donde u cambia de signo
        return None, _evaluar_por_tramos(f, an, bn)
    F = antiderivada(f)
    return F, sp.simplify(F.subs(x, bn) - F.subs(x, an))


def ceros_en(f, a, b):
    """Ceros de f en (a, b): exactos si sympy los encuentra; si no, por cambio de signo."""
    sol = sp.solveset(f, x, sp.Interval.open(sp.nsimplify(a), sp.nsimplify(b)))
    if isinstance(sol, sp.FiniteSet):
        return sorted(sol, key=float)
    g = numerica(f); xs = np.linspace(float(a), float(b), 2001); y = g(xs)
    return [sp.Float(brentq(lambda s: float(g(np.array([s]))[0]), p, q), 12)
            for p, q, yp, yq in zip(xs, xs[1:], y, y[1:]) if yp * yq < 0]


def _evaluar_por_tramos(f, a, b):
    """∫_a^b de f con valores absolutos: en cada tramo, |u| se sustituye por u o por -u."""
    lo, hi, signo = (a, b, 1) if a < b else (b, a, -1)
    cortes = {lo, hi}
    for u in [e.args[0] for e in f.atoms(sp.Abs)]:
        cortes |= set(ceros_en(u, lo, hi))
    cortes = sorted(cortes, key=float)
    total = 0
    for p, q in zip(cortes, cortes[1:]):
        m = (p + q) / 2
        g = f.replace(sp.Abs, lambda u: u if float(u.subs(x, m)) >= 0 else -u)
        G = antiderivada(g)
        total += G.subs(x, q) - G.subs(x, p)
    return sp.simplify(signo * total)


def punto_medio(f, a, b, m=1000):
    xs = np.linspace(a, b, m + 1)
    return float(numerica(f)((xs[:-1] + xs[1:]) / 2).sum() * (b - a) / m)


def calculadora_evaluacion(f_txt, a, b, n_cifras):
    try:
        f = parsear(f_txt)
        F, valor = evaluar(f, a, b)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"F(x) = {F}" if F is not None else "f tiene valores absolutos: se integró por tramos")
    print(f"[F(x)] de {a} a {b} = {valor} ≈ {cifras(valor, n_cifras)}   ({n_cifras} cifras significativas, redondeado)")
    print(f"comprobación con M_1000 ≈ {cifras(punto_medio(f, a, b), n_cifras)}")


widgets.interact(calculadora_evaluacion,
    f_txt=widgets.Text(value="x**2", description="f(x) ="),
    a=widgets.FloatText(value=0, description="a"), b=widgets.FloatText(value=3, description="b"),
    n_cifras=widgets.IntSlider(value=6, min=1, max=10, description="cifras sig."));

# %%
# Casos de prueba de la sección 2 (resultado conocido)
PRUEBAS_2 = [
    ("∫_0^π sen x dx = 2", lambda: evaluar(sp.sin(x), 0, sp.pi)[1] == 2),
    ("∫_1^e dx/x = 1", lambda: evaluar(1 / x, 1, sp.E)[1] == 1),
    ("∫_0^1 e^(2x) dx = (e^2 - 1)/2", lambda: sp.simplify(evaluar(sp.exp(2 * x), 0, 1)[1] - (sp.E**2 - 1) / 2) == 0),
    ("1/x^2 en [-1, 1] se rechaza (no continua)", lambda: _rechaza(lambda: evaluar(1 / x**2, -1, 1))),
    ("e^(-x^2) en [0, 1]: sin antiderivada elemental, se avisa", lambda: _rechaza(lambda: evaluar(sp.exp(-x**2), 0, 1))),
    ("|x - 1| en [0, 3]: 5/2, por tramos", lambda: evaluar(sp.Abs(x - 1), 0, 3)[1] == sp.Rational(5, 2)),
]
for nombre, prueba in PRUEBAS_2:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 2).** Calcula a mano $\int_1^4(2x-1)\,dx$.

# %%
mi_valor = None     # escribe un número

if mi_valor is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    F, ref = evaluar(2 * x - 1, 1, 4)
    print("La calculadora da: F(x) =", F, "  valor =", ref)
    print("coinciden" if cerca(mi_valor, ref, 1e-9) else "NO coinciden: F(x) = x^2 - x; evalúa en 4 y en 1 y resta en ese orden")
