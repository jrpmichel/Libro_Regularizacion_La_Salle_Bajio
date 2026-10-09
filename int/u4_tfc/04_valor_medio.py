# ID: INT-U4-NB04
# Notebook: int/u4_tfc.ipynb · sección 4.4 teorema del valor medio para integrales
# Repositorio: int/u4_tfc/04_valor_medio.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 4. Valor medio de una función
#
# $$\bar f=\frac{1}{b-a}\int_a^bf(x)\,dx,\qquad f(c)=\bar f\ \text{para algún } c\in[a,b]\ \text{si } f \text{ es continua}.$$
#
# La calculadora da $\bar f$, busca **todos** los $c$ de $[a,b]$ con $f(c)=\bar f$ y dibuja el rectángulo de altura $\bar f$. También muestra el promedio de los extremos, $\frac{f(a)+f(b)}{2}$, que coincide con $\bar f$ siempre que $f$ es lineal y, en otros casos, solo por casualidad.

# %%
def valor_medio(f, a, b):
    """(valor medio exacto, lista de c en [a, b] o 'todo el intervalo')."""
    if not a < b:
        raise ValueError("el intervalo debe cumplir a < b.")
    fbar = sp.simplify(evaluar(f, a, b)[1] / (sp.nsimplify(b) - sp.nsimplify(a)))
    escala = max(1.0, float(np.max(np.abs(numerica(f)(np.linspace(float(a), float(b), 401))))))
    if abs(float(fbar)) < 1e-12 * escala:      # ruido de redondeo: el valor medio es cero
        fbar = sp.Integer(0)
    if sp.simplify(f - fbar) == 0:
        return fbar, "todo el intervalo"
    g = numerica(f - fbar)
    xs = np.linspace(float(a), float(b), 4001); y = g(xs)
    cs = [brentq(lambda s: float(g(np.array([s]))[0]), p, q) for p, q, yp, yq in zip(xs, xs[1:], y, y[1:]) if yp * yq < 0]
    cs += [float(p) for p, yp in zip(xs, y) if abs(yp) < 1e-12 * escala]     # incluye los extremos
    return fbar, sorted(set(round(c, 9) for c in cs))


def calculadora_valor_medio(f_txt, a, b, n_cifras):
    try:
        f = parsear(f_txt)
        fbar, cs = valor_medio(f, a, b)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"valor medio = {fbar} ≈ {cifras(fbar, n_cifras)}   ({n_cifras} cifras significativas, redondeado)")
    print("puntos c:", cs if isinstance(cs, str) else ", ".join(cifras(c, n_cifras) for c in cs))
    extremos = (valor_real(f, a) + valor_real(f, b)) / 2
    print(f"promedio de los extremos ≈ {cifras(0 if abs(extremos) < 1e-12 else extremos, n_cifras)}")
    xs = np.linspace(a, b, 400)
    fig, ax = plt.subplots(figsize=(5.5, 3))
    ax.fill_between([a, b], 0, [float(fbar)] * 2, color="0.85", label="rectángulo de altura f̄")
    ax.plot(xs, numerica(f)(xs), color="navy", label="f(x)")
    ax.axhline(float(fbar), color="red", ls="--")
    if not isinstance(cs, str):
        ax.plot(cs, [float(fbar)] * len(cs), "o", color="red")
    ax.legend(); plt.show()


widgets.interact(calculadora_valor_medio,
    f_txt=widgets.Text(value="2 + sin(x)", description="f(x) ="),
    a=widgets.FloatText(value=0, description="a"), b=widgets.FloatText(value=4, description="b"),
    n_cifras=widgets.IntSlider(value=5, min=1, max=10, description="cifras sig."));

# %%
# Casos de prueba de la sección 4 (resultado conocido)
PRUEBAS_4 = [
    ("x^2 en [0, 3]: valor medio 3, c = √3",
     lambda: valor_medio(x**2, 0, 3)[0] == 3 and cerca(valor_medio(x**2, 0, 3)[1][0], math.sqrt(3), 1e-9)),
    ("sen x en [0, π]: 2/π y dos puntos c",
     lambda: valor_medio(sp.sin(x), 0, sp.pi)[0] == 2 / sp.pi and len(valor_medio(sp.sin(x), 0, sp.pi)[1]) == 2),
    ("3x en [0, 2]: valor medio 3 en c = 1", lambda: valor_medio(3 * x, 0, 2)[0] == 3 and cerca(valor_medio(3 * x, 0, 2)[1][0], 1)),
    ("constante 5 en [1, 4]: todo el intervalo", lambda: valor_medio(5 + 0 * x, 1, 4) == (5, "todo el intervalo")),
    ("2 + sen x en [0, 4]: c ≈ 0.4262 y 2.7154",
     lambda: [round(c, 4) for c in valor_medio(2 + sp.sin(x), 0, 4)[1]] == [0.4262, 2.7154]),
]
for nombre, prueba in PRUEBAS_4:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 4).** Calcula a mano el valor medio de $f(x)=2x+1$ en $[0,4]$ y el punto $c$.

# %%
mi_media = None   # escribe un número
mi_c = None       # escribe un número

if mi_media is None or mi_c is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    fbar, cs = valor_medio(2 * x + 1, 0, 4)
    print("La calculadora da: valor medio =", fbar, "  c =", cs)
    print("coinciden" if cerca(mi_media, fbar) and cerca(mi_c, cs[0]) else
          "NO coinciden: la integral vale 20 y el intervalo mide 4; con f lineal, c es el punto medio")
