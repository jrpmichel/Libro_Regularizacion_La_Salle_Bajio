# ID: INT-U6-NB03
# Notebook: int/u6_aplicaciones.ipynb · sección 6.3 longitud de arco
# Repositorio: int/u6_aplicaciones/03_arco.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6.3 Longitud de arco
#
# $$L=\int_a^b\sqrt{1+\big(f'(x)\big)^2}\,dx,\qquad L_n=\sum_{k=1}^{n}\sqrt{(\Delta x)^2+(\Delta y_k)^2}\quad(\text{poligonal de }n\text{ segmentos iguales en }x).$$
#
# Cada segmento de la poligonal es una cuerda, más corta que el arco que une sus extremos, así que $L_n\le L$ y $L_n\to L$ cuando $n$ crece. La integral casi nunca tiene antiderivada elemental; cuando no la tiene, se calcula numéricamente.
#
# El ejemplo exacto usa $y=x^{3/2}$ en $[0,4]$, cuyo integrando $\sqrt{1+\tfrac94x}$ sí tiene antiderivada.

# %%
y_ej = x**sp.Rational(3, 2)
integrando_ej = sp.sqrt(1 + sp.diff(y_ej, x)**2)
L_ej = sp.integrate(integrando_ej, (x, 0, 4))
print("√(1 + f'²) =", integrando_ej, "   L =", L_ej, "≈", cifras(L_ej, 6))

# %% [markdown]
# **Calculadora de longitud de arco.** Escribe $f$, los extremos y el número $n$ de segmentos. La calculadora intenta la integral exacta con `sympy` (con límite de tiempo); si no sale, la calcula con `quad` y te dice cuál usó. Mueve $n$ para ver cómo la poligonal se acerca a la curva por debajo.

# %%
def longitud_arco(f, a, b):
    """L = ∫_a^b √(1 + f'²) dx. Devuelve un dict con L, el integrando y el método (exacto con sympy o quad)."""
    f = sp.sympify(f)
    a, b = exacto(a, "a"), exacto(b, "b")
    if not a < b:
        raise ValueError("a debe ser menor que b.")
    revisa_intervalo(f, a, b, "f")
    df = sp.diff(f, x)
    if f.has(sp.Heaviside) or df.has(sp.DiracDelta):
        raise ValueError("f tiene un salto (un escalón): la longitud de arco pide una curva continua con f' continua, "
                         "al menos por tramos. Mide cada tramo continuo por separado.")
    try:                                   # factor: a veces 1 + f'² es un cuadrado perfecto (catenaria)
        integrando = con_limite(lambda: sp.sqrt(sp.factor(sp.together(1 + df**2))), 5)
    except ValueError:
        integrando = sp.sqrt(1 + df**2)
    L, como = integra(integrando, a, b)
    return {"L": L, "integrando": integrando, "metodo": como}


def poligonal(f, a, b, n):
    """Longitud de la poligonal de n segmentos con vértices sobre y = f(x), igualmente espaciados en x."""
    if int(n) != n or n < 1:
        raise ValueError("n debe ser un entero mayor o igual que 1.")
    a, b = float(exacto(a, "a")), float(exacto(b, "b"))
    if not a < b:
        raise ValueError("a debe ser menor que b.")
    xs = np.linspace(a, b, int(n) + 1)
    ys = numerica(sp.sympify(f))(xs)
    if not np.all(np.isfinite(ys)):
        raise ValueError("f no tiene valor real en algún vértice de la poligonal: revisa el dominio.")
    return float(np.sum(np.hypot(np.diff(xs), np.diff(ys))))


def calculadora_arco(f_txt, a_txt, b_txt, n):
    try:
        f = parsear(f_txt)
        a, b = exacto(a_txt, "a"), exacto(b_txt, "b")
        r = longitud_arco(f, a, b)
        Ln = poligonal(f, a, b, n)
    except ERRORES as err:
        print("Revisa la entrada:", explica(err)); return
    print(f"({CIFRAS} cifras significativas, redondeadas; si x y y están en metros, L sale en m)")
    print(f"integrando √(1 + f'²) = {r['integrando']}")
    print(f"{'L (integral)':<20} = {muestra(r['L'])} u   ({COMO[r['metodo']]})")
    print(f"{f'L_{n} (poligonal)':<20} = {cifras(Ln, CIFRAS)} u")
    print(f"{f'diferencia L - L_{n}':<20} = {cifras(float(r['L']) - Ln, 3)} u   (3 cifras)")
    try:
        xs = np.linspace(float(a), float(b), 600)
        xp = np.linspace(float(a), float(b), n + 1)
        ys, yp = para_graficar(f, xs), para_graficar(f, xp)
        fig, ax = plt.subplots(figsize=(5.5, 3.2))
        ax.plot(xs, ys, color="black", ls="-")
        ax.plot(xp, yp, color="0.35", ls="--", marker="o" if n <= 32 else None, ms=3)
        rotula(ax, xs[-1], ys[-1], "curva")
        k = max(0, n // 2 - 1)
        rotula(ax, (xp[k] + xp[k + 1]) / 2, (yp[k] + yp[k + 1]) / 2, f"poligonal, n = {n}", xytext=(6, -10), va="top")
        ax.set_xlabel("x"); ax.set_title("arco (continua) y poligonal inscrita (discontinua)", fontsize=9); plt.show()
    except ERRORES:
        print("(no pude dibujar la gráfica con estos datos)")


widgets.interact(calculadora_arco,
    f_txt=widgets.Text(value="x^(3/2)", description="f(x) =", continuous_update=False),
    a_txt=widgets.Text(value="0", description="a =", continuous_update=False),
    b_txt=widgets.Text(value="4", description="b =", continuous_update=False),
    n=widgets.IntSlider(value=4, min=1, max=64, description="n", continuous_update=False));

# %%
# Casos de prueba de la sección 6.3 (resultado conocido)
PRUEBAS_3 = [
    ("y = x^(3/2) en [0, 4]: (80√10 - 8)/27 ≈ 9.07342, exacta",
     lambda: (lambda r: r["metodo"] == "exacta" and sp.simplify(r["L"] - (80 * sp.sqrt(10) - 8) / 27) == 0
              and cifras(r["L"], 6) == "9.07342")(longitud_arco(x**sp.Rational(3, 2), 0, 4))),
    ("cable parabólico y = 4·20·x²/200² en [-100, 100]: 205.212 m",
     lambda: round(float(longitud_arco(4 * 20 * x**2 / 200**2, -100, 100)["L"]), 3) == 205.212),
    ("catenaria y = (e^x + e^(-x))/2 en [-1, 1]: e - 1/e",
     lambda: cerca(longitud_arco((sp.exp(x) + sp.exp(-x)) / 2, -1, 1)["L"], math.e - 1 / math.e, 1e-12)),
    ("poligonal de y = x^(3/2) con n = 64: a menos de 1e-3 de la exacta, y por debajo",
     lambda: 0 <= float((80 * sp.sqrt(10) - 8) / 27) - poligonal(x**sp.Rational(3, 2), 0, 4, 64) < 1e-3),
    ("poligonal de una recta con n = 1: exacta (y = 3x + 1 en [0, 4] → 4√10)",
     lambda: cerca(poligonal(3 * x + 1, 0, 4, 1), 4 * math.sqrt(10), 1e-12)),
    ("y = sen x en [0, π]: sin antiderivada elemental, se usa quad (3.82020)",
     lambda: (lambda r: r["metodo"] == "numérica" and cifras(r["L"], 6) == "3.82020")(longitud_arco(sp.sin(x), 0, sp.pi))),
    ("n = 0 se rechaza", lambda: _rechaza(lambda: poligonal(x**2, 0, 1, 0))),
]
for nombre, prueba in PRUEBAS_3:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 6.3).** Calcula a mano la longitud de la recta $y=2x$ entre $x=0$ y $x=3$ y escribe el número.

# %%
mi_valor = None     # escribe un número (también sirve "3*sqrt(5)")

if mi_valor is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    try:
        mio = a_numero(mi_valor)
    except ValueError as err:
        print("Revisa tu número:", err)
    else:
        ref = longitud_arco(2 * x, 0, 3)["L"]
        print("La calculadora da:", muestra(ref), "u")
        print("coinciden a 3 cifras" if coincide(mio, ref) else
              "NO coinciden: con f' = 2 el integrando es √(1 + 4) = √5, constante; multiplica por la longitud del intervalo")
