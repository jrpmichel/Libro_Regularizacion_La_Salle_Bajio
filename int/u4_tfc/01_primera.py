# ID: INT-U4-NB01
# Notebook: int/u4_tfc.ipynb · sección 4.1 primera parte: la derivada de la integral
# Repositorio: int/u4_tfc/01_primera.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Primera parte: la derivada de la integral
#
# Si $f$ es continua en $[a,b]$ y $G(x)=\int_a^x f(t)\,dt$, entonces $G'(x)=f(x)$. Con límites que dependen de $x$:
#
# $$\frac{d}{dx}\int_{u(x)}^{v(x)}f(t)\,dt=f\big(v(x)\big)\,v'(x)-f\big(u(x)\big)\,u'(x).$$
#
# La primera calculadora construye $G$ con trapecios y compara su derivada numérica con $f$. La segunda deriva una integral con límites variables y muestra cada término de la regla de la cadena. Escribe $f$ con la variable `x`; la calculadora la lee como $f(t)$.

# %%
def acumulada(f, a, b, m=2000):
    """(t, G) con G(t) = ∫_a^t f, acumulada con trapecios en m franjas."""
    if not a < b:
        raise ValueError("el intervalo debe cumplir a < b.")
    if not continua_en(f, a, b):
        raise ValueError("f no es continua en todo [a, b]; la primera parte del teorema pide continuidad.")
    F = numerica(f)
    t = np.linspace(a, b, m + 1)
    y = F(t)
    G = np.concatenate([[0.0], np.cumsum((y[:-1] + y[1:]) / 2 * np.diff(t))])
    return t, G


def derivada_integral(f, u, v):
    """d/dx de ∫_{u(x)}^{v(x)} f(t) dt y sus dos términos."""
    alto = f.subs(x, v) * sp.diff(v, x)
    bajo = f.subs(x, u) * sp.diff(u, x)
    return sp.simplify(alto - bajo), sp.simplify(alto), sp.simplify(bajo)


def calculadora_acumulada(f_txt, a, b, n_cifras):
    try:
        f = parsear(f_txt)
        t, G = acumulada(f, a, b)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    dG = np.gradient(G, t[1] - t[0])
    print(f"G({b}) = ∫ de {a} a {b} ≈ {cifras(G[-1], n_cifras)}   ({n_cifras} cifras significativas, redondeado)")
    print("máx |G' - f| en puntos interiores ≈", cifras(np.max(np.abs(dG[1:-1] - numerica(f)(t[1:-1]))), 2))
    fig, ax = plt.subplots(figsize=(5.5, 3))
    ax.plot(t, numerica(f)(t), color="navy", label="f(t)")
    ax.plot(t, G, color="red", label="G(x)")
    ax.axhline(0, color="black", lw=0.6); ax.set_xlabel("x"); ax.legend(); plt.show()


def calculadora_limites(f_txt, u_txt, v_txt):
    try:
        f, u, v = parsear(f_txt), parsear(u_txt), parsear(v_txt)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    total, alto, bajo = derivada_integral(f, u, v)
    if sp.calculus.util.continuous_domain(f, x, sp.S.Reals) != sp.S.Reals:
        print("Aviso: f no es continua en toda la recta; la fórmula vale solo si f es continua entre los dos límites.")
    print(f"d/dx ∫ de {u} a {v} de f(t) dt, con f(t) = {f.subs(x, sp.Symbol('t'))}")
    print("  término del límite superior: f(v)·v' =", alto)
    print("  término del límite inferior: f(u)·u' =", bajo)
    print("  derivada =", total)


widgets.interact(calculadora_acumulada,
    f_txt=widgets.Text(value="2 + cos(x)", description="f ="),
    a=widgets.FloatText(value=0, description="a"), b=widgets.FloatText(value=4, description="b"),
    n_cifras=widgets.IntSlider(value=5, min=1, max=8, description="cifras sig."));
widgets.interact(calculadora_limites,
    f_txt=widgets.Text(value="cos(x)", description="f ="),
    u_txt=widgets.Text(value="0", description="u(x) ="), v_txt=widgets.Text(value="x**2", description="v(x) ="));

# %%
# Casos de prueba de la sección 1 (resultado conocido)
PRUEBAS_1 = [
    ("d/dx ∫_0^x t^2 dt = x^2", lambda: misma_funcion(derivada_integral(x**2, 0 * x, x)[0], x**2)),
    ("d/dx ∫_1^{x^2} dt/t = 2/x", lambda: misma_funcion(derivada_integral(1 / x, 1 + 0 * x, x**2)[0], 2 / x)),
    ("d/dx ∫_x^5 t dt = -x", lambda: misma_funcion(derivada_integral(x, x, 5 + 0 * x)[0], -x)),
    ("G de cos en [0, π]: G(π/2) ≈ 1", lambda: cerca(np.interp(math.pi / 2, *acumulada(sp.cos(x), 0, math.pi)), 1, 1e-6)),
    ("e^(-t^2) en [0, 2]: máx |G' - f| < 1e-4",
     lambda: (lambda t, G: np.max(np.abs(np.gradient(G, t[1] - t[0])[1:-1] - np.exp(-t[1:-1]**2))) < 1e-4)(
         *acumulada(sp.exp(-x**2), 0, 2))),
]
for nombre, prueba in PRUEBAS_1:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 1).** Deriva a mano $H(x)=\int_0^{3x}e^{t^2}\,dt$.

# %%
mi_H = None     # escribe H'(x) como texto, con x como variable

if mi_H is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    ref = derivada_integral(sp.exp(x**2), 0 * x, 3 * x)[0]
    print("La calculadora da:", ref)
    try:
        ok = misma_funcion(parsear(str(mi_H)), ref)
    except ValueError as err:
        print(err); ok = False
    print("coinciden" if ok else "NO coinciden: evalúa e^(t^2) en el límite superior, 3x, y multiplica por la derivada de 3x")
