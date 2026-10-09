# ID: INT-U1-NB02
# Notebook: int/u1_area.ipynb · sección 1.2 sumas inferior y superior
# Repositorio: int/u1_area/02_rectangulos.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Aproximación con rectángulos: sumas inferior y superior
#
# Partición regular: $\Delta x=\dfrac{b-a}{n}$, $x_k=a+k\,\Delta x$. En cada subintervalo $[x_{k-1},x_k]$, $m_k$ es el **mínimo** de $f$ y $M_k$ el **máximo**:
#
# $$s_n=\sum m_k\,\Delta x\ \ (\text{inscritos})\qquad\le\qquad A\qquad\le\qquad S_n=\sum M_k\,\Delta x\ \ (\text{circunscritos}).$$
#
# Si $f$ es **monótona**, el mínimo y el máximo están en los extremos de cada subintervalo y $S_n-s_n=|f(b)-f(a)|\,\Delta x$. Si **no** lo es, pueden estar dentro: la calculadora busca también los puntos donde $f'=0$.

# %%
def _valor(f, c):
    v = complex(sp.N(f.subs(x, c)))
    if abs(v.imag) > 1e-12 or not math.isfinite(v.real):
        raise ValueError(f"f no tiene valor real en x = {cifras(float(c), 6)}: revisa el dominio.")
    return v.real


def extremos_en(f, a, b):
    """(mínimo, máximo) de f en [a, b]: extremos del subintervalo y puntos interiores con f' = 0."""
    cand = [sp.nsimplify(a), sp.nsimplify(b)]
    crit = sp.solveset(sp.diff(f, x), x, domain=sp.Interval.open(sp.nsimplify(a), sp.nsimplify(b)))
    if isinstance(crit, sp.FiniteSet):
        cand += list(crit)
    elif crit is not sp.S.EmptySet:              # sympy no pudo listarlos: muestreo fino como respaldo
        cand += list(np.linspace(float(a), float(b), 401))
    vals = [_valor(f, c) for c in cand]
    return min(vals), max(vals)


def sumas_inf_sup(f, a, b, n):
    """Suma inferior, suma superior y lista de (m_k, M_k) con partición regular de n subintervalos."""
    if not a < b:
        raise ValueError("el intervalo debe cumplir a < b.")
    intervalo = sp.Interval(sp.nsimplify(a), sp.nsimplify(b))
    if not intervalo.is_subset(sp.calculus.util.continuous_domain(f, x, intervalo)):
        raise ValueError("f no es continua en todo [a, b] (tiene un salto, una asíntota o sale de su dominio); "
                         "las sumas inferior y superior de esta sección suponen f continua.")
    if int(n) != n or n < 1:
        raise ValueError("n debe ser un entero positivo.")
    n = int(n)
    xs = [sp.nsimplify(a) + k * (sp.nsimplify(b) - sp.nsimplify(a)) / n for k in range(n + 1)]
    dx = float(xs[1] - xs[0])
    mM = [extremos_en(f, xs[k], xs[k + 1]) for k in range(n)]
    if min(m for m, _ in mM) < -1e-12:
        raise ValueError("f toma valores negativos en [a, b]; esta sección trabaja con regiones sobre el eje.")
    return sum(m for m, _ in mM) * dx, sum(M for _, M in mM) * dx, mM


def monotona(f, a, b):
    """True si f' no cambia de signo en (a, b) (revisión en 400 puntos)."""
    d = sp.lambdify(x, sp.diff(f, x), "numpy")
    v = np.asarray(d(np.linspace(a, b, 402)[1:-1]), dtype=float) * np.ones(400)
    return bool(np.all(v >= -1e-12) or np.all(v <= 1e-12))


def n_para_tolerancia(f, a, b, tol):
    """Rectángulos que garantizan S_n - s_n <= tol si f es monótona (propiedad 1.6b)."""
    if not monotona(f, a, b):
        raise ValueError("f no es monótona en [a, b]: la fórmula |f(b) - f(a)|(b - a)/n no aplica.")
    salto = abs(_valor(f, b) - _valor(f, a)) * (b - a)
    return max(1, math.ceil(salto / tol - 1e-12))


def calculadora_rectangulos(f_txt, a, b, n, n_cifras, tolerancia):
    try:
        f = parsear(f_txt)
        s, S, mM = sumas_inf_sup(f, a, b, n)
    except (ValueError, TypeError) as err:
        print("Revisa la entrada:", err); return
    print(f"s_{n} ≈ {cifras(s, n_cifras)}   S_{n} ≈ {cifras(S, n_cifras)}   diferencia ≈ {cifras(S - s, n_cifras)}"
          f"   promedio ≈ {cifras((s + S) / 2, n_cifras)}   ({n_cifras} cifras significativas, redondeado)")
    try:
        print(f"para S_n - s_n <= {tolerancia}: n >= {n_para_tolerancia(f, a, b, tolerancia)}")
    except ValueError as err:
        print("tolerancia:", err)
    fx = sp.lambdify(x, f, "numpy")
    xx = np.linspace(a, b, 400)
    fig, ax = plt.subplots(figsize=(5.5, 3))
    dx = (b - a) / n
    for k, (m, M) in enumerate(mM):
        ax.add_patch(plt.Rectangle((a + k * dx, 0), dx, M, fill=False, hatch="///", edgecolor="gray"))
        ax.add_patch(plt.Rectangle((a + k * dx, 0), dx, m, color="0.8"))
    ax.plot(xx, np.asarray(fx(xx), dtype=float) * np.ones_like(xx), color="navy")
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xlabel("x"); ax.set_ylabel("f(x)")
    ax.set_title(f"inscritos (gris) y circunscritos (achurado), n = {n}"); plt.show()


widgets.interact(calculadora_rectangulos,
    f_txt=widgets.Text(value="x**2", description="f(x) ="),
    a=widgets.FloatText(value=0, description="a"), b=widgets.FloatText(value=1, description="b"),
    n=widgets.IntSlider(value=4, min=1, max=60, description="n"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."),
    tolerancia=widgets.FloatText(value=0.01, description="tolerancia"));

# %%
# Casos de prueba de la sección 2 (resultado conocido)
_s, _S, _ = sumas_inf_sup(x**2, 0, 1, 4)
_si, _Si, _ = sumas_inf_sup(sp.sin(x), 0, sp.pi, 3)
PRUEBAS_2 = [
    ("x^2 en [0, 1], n = 4: s = 7/32 y S = 15/32", lambda: cerca(_s, 7 / 32) and cerca(_S, 15 / 32)),
    ("1/x en [1, 2], n = 4: s ≈ 0.63452, S ≈ 0.75952",
     lambda: cerca(sumas_inf_sup(1 / x, 1, 2, 4)[0], 0.634523809, 1e-8) and cerca(sumas_inf_sup(1 / x, 1, 2, 4)[1], 0.759523809, 1e-8)),
    ("constante 3 en [0, 2], n = 5: s = S = 6", lambda: cerca(sumas_inf_sup(3 + 0 * x, 0, 2, 5)[0], 6)
                                                    and cerca(sumas_inf_sup(3 + 0 * x, 0, 2, 5)[1], 6)),
    ("sin x en [0, π], n = 3: máximo interior; s ≈ 0.9069, S ≈ 2.861",
     lambda: cerca(_si, 0.906899682, 1e-8) and cerca(_Si, 2.860996915, 1e-8)),
    ("x^2 en [0, 1]: S_n - s_n <= 0.001 exige n >= 1000", lambda: n_para_tolerancia(x**2, 0, 1, 0.001) == 1000),
    ("a >= b, f negativa o f discontinua se rechazan", lambda: _rechaza(lambda: sumas_inf_sup(x, 2, 1, 4))
                                                and _rechaza(lambda: sumas_inf_sup(x - 2, 0, 1, 4))
                                                and _rechaza(lambda: sumas_inf_sup(1 / (x - sp.Rational(1, 2))**2, 0, 1, 3))),
]
for nombre, prueba in PRUEBAS_2:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 2).** Calcula a mano $s_4$ y $S_4$ para $f(x)=x$ en $[0,2]$. Comprueba que el área exacta, que conoces por geometría, queda entre tus dos sumas.

# %%
mi_s4 = None     # escribe un número
mi_S4 = None     # escribe un número

if mi_s4 is None or mi_S4 is None:      # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    s4, S4, _ = sumas_inf_sup(x, 0, 2, 4)
    print("La calculadora da: s_4 =", cifras(s4, 4), " S_4 =", cifras(S4, 4), " (área exacta: triángulo de base 2 y altura 2)")
    print("coinciden" if cerca(mi_s4, s4, 1e-6) and cerca(mi_S4, S4, 1e-6) else
          "NO coinciden: con Δx = 0.5, la suma inferior usa 0, 0.5, 1 y 1.5 como alturas")
