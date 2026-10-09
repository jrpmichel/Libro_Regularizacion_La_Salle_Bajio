# ID: INT-U6-NB04
# Notebook: int/u6_aplicaciones.ipynb · sección 6.4 centroides y momentos de inercia de área
# Repositorio: int/u6_aplicaciones/04_centroides.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6.4 Centroides y momentos de inercia de área
#
# Para la región entre $g$ (abajo) y $f$ (arriba) en $[a,b]$:
#
# $$A=\int_a^b(f-g)\,dx,\qquad \bar x=\frac1A\int_a^b x\,(f-g)\,dx,\qquad \bar y=\frac1A\int_a^b\frac{f^2-g^2}{2}\,dx,$$
#
# $$I_x=\int_a^b\frac{f^3-g^3}{3}\,dx\ \ (\text{respecto al eje }x),\qquad I_{\bar x}=I_x-A\,\bar y^{\,2}\ \ (\text{teorema de Steiner}).$$
#
# Para una sección hecha de rectángulos $b_i\times h_i$ con centro a la altura $y_i$: $\bar y=\dfrac{\sum A_iy_i}{\sum A_i}$ e $I_{\bar x}=\sum\left(\dfrac{b_ih_i^3}{12}+A_i\,d_i^2\right)$, con $d_i=y_i-\bar y$. El término $A_id_i^2$ traslada la inercia de cada rectángulo al eje centroidal de la sección; sin él la viga parece mucho más flexible de lo que es.
#
# El ejemplo exacto usa la región bajo $y=4-x^2$ en $[0,2]$.

# %%
f_ej = 4 - x**2
A_ej = sp.integrate(f_ej, (x, 0, 2))
xb_ej = sp.integrate(x * f_ej, (x, 0, 2)) / A_ej
yb_ej = sp.integrate(f_ej**2 / 2, (x, 0, 2)) / A_ej
print("A =", A_ej, "  x̄ =", xb_ej, "  ȳ =", yb_ej)

# %% [markdown]
# **Calculadora de una región entre curvas.** Escribe $f$ (arriba), $g$ (abajo) y los extremos, con longitudes en cm. La calculadora da el área, el centroide, $I_x$ respecto al eje $x$ e $I_{\bar x}$ respecto al eje horizontal que pasa por el centroide, y los dibuja.

# %%
def region(f, g, a, b):
    """Área, centroide (x̄, ȳ), I_x respecto al eje x e I_c respecto al eje centroidal horizontal (Steiner)
    de la región entre g (abajo) y f (arriba) en [a, b]. Devuelve un dict. Si una región no acotada tiene área
    finita pero x̄, ȳ o los momentos divergen (como bajo 1/√x en [0, 1]), esos valores son None y la clave
    "divergen" dice cuáles."""
    f, g = sp.sympify(f), sp.sympify(g)
    a, b = exacto(a, "a"), exacto(b, "b")
    if not a < b:
        raise ValueError("a debe ser menor que b.")
    revisa_intervalo(f, a, b, "f"); revisa_intervalo(g, a, b, "g")
    xs = malla_interior(a, b)
    dif = numerica(f)(xs) - numerica(g)(xs)
    if np.any(dif < -1e-9 * max(1.0, float(np.max(np.abs(dif))))):
        raise ValueError(f"en x ≈ {cifras(xs[int(np.argmin(dif))], 4)} la curva g queda arriba de f. Escribe arriba "
                         "la curva superior (f) o separa la región en el cruce.")
    integ = Integrador()
    A = integ(sp.expand(f - g), a, b)
    if float(A) <= 0:
        raise ValueError("la región tiene área cero: f y g coinciden en todo el intervalo.")
    def _opcional(h):
        try:
            return integ(sp.expand(h), a, b)
        except ValueError:                          # diverge: el resto del resultado sigue valiendo
            return None
    Mx, My, Ix = _opcional(x * (f - g)), _opcional((f**2 - g**2) / 2), _opcional((f**3 - g**3) / 3)
    xb = simplifica(Mx / A) if Mx is not None else None
    yb = simplifica(My / A) if My is not None else None
    Ic = simplifica(Ix - A * yb**2) if Ix is not None and yb is not None else None
    divergen = [n for n, v in (("x̄", xb), ("ȳ", yb), ("I_x", Ix), ("I_c", Ic)) if v is None]
    return {"A": A, "x_bar": xb, "y_bar": yb, "I_x": Ix, "I_c": Ic, "metodo": integ.metodo, "divergen": divergen}


def calculadora_region(f_txt, g_txt, a_txt, b_txt):
    try:
        f, g = parsear(f_txt), parsear(g_txt)
        a, b = exacto(a_txt, "a"), exacto(b_txt, "b")
        r = region(f, g, a, b)
    except ERRORES as err:
        print("Revisa la entrada:", explica(err)); return
    def _v(clave, unidad):
        return f"{muestra(r[clave])} {unidad}" if r[clave] is not None else "no existe: su integral diverge"
    print(f"({CIFRAS} cifras significativas, redondeadas; longitudes en cm; {COMO[r['metodo']]})")
    print(f"A   = {muestra(r['A'])} cm²")
    print(f"x̄   = {_v('x_bar', 'cm')}     ȳ = {_v('y_bar', 'cm')}")
    print(f"I_x = {_v('I_x', 'cm⁴')}   (respecto al eje x)")
    print(f"I_c = {_v('I_c', 'cm⁴')}   (respecto al eje centroidal horizontal, I_x - A ȳ²)")
    if r["divergen"]:
        print("La región no está acotada cerca de un extremo: su área es finita, pero", ", ".join(r["divergen"]),
              "no, porque usan f² o f³, que crecen más rápido cerca de la asíntota.")
    try:
        xs = np.linspace(float(a), float(b), 400)
        yf, yg = para_graficar(f, xs), para_graficar(g, xs)
        fig, ax = plt.subplots(figsize=(5.5, 3.2))
        ax.fill_between(xs, yg, yf, color="0.85")
        ax.plot(xs, yf, color="black", ls="-"); ax.plot(xs, yg, color="black", ls="--")
        k1, k2 = len(xs) // 4, len(xs) * 3 // 4
        rotula(ax, xs[k1], yf[k1], "f", xytext=(0, 9), ha="center")
        rotula(ax, xs[k2], yg[k2], "g", xytext=(0, -9), ha="center")
        if r["x_bar"] is not None and r["y_bar"] is not None:
            xb, yb = float(r["x_bar"]), float(r["y_bar"])
            ax.axhline(yb, color="0.3", ls="-.", lw=0.8); rotula(ax, xs[-1], yb, "eje centroidal", xytext=(6, 6))
            ax.plot([xb], [yb], marker="+", color="black", ms=14, mew=2)
            rotula(ax, xb, yb, f"C ({cifras(xb, 3)}, {cifras(yb, 3)})", xytext=(8, -10))
        ax.set_xlabel("x (cm)"); ax.set_ylabel("y (cm)"); ax.set_aspect("equal", adjustable="datalim")
        ax.set_title("C: centroide (coordenadas con 3 cifras)", fontsize=9); plt.show()
    except ERRORES:
        print("(no pude dibujar la gráfica con estos datos)")


widgets.interact(calculadora_region,
    f_txt=widgets.Text(value="4 - x^2", description="f(x) =", continuous_update=False),
    g_txt=widgets.Text(value="0", description="g(x) =", continuous_update=False),
    a_txt=widgets.Text(value="0", description="a =", continuous_update=False),
    b_txt=widgets.Text(value="2", description="b =", continuous_update=False));

# %% [markdown]
# **Calculadora de una sección compuesta de rectángulos.** Escribe cada rectángulo como `b, h, y_centro` en cm (base, altura y altura de su centro); deja vacío el que no uses. Los valores de ejemplo forman una sección T: patín de $12\times2$ con centro en $y=11$ y alma de $2\times10$ con centro en $y=5$. La calculadora suma los rectángulos con el teorema de Steiner y también muestra la suma sin el término $A_id_i^2$, para que veas el tamaño del error.

# %%
def seccion_compuesta(rectangulos):
    """Sección hecha de rectángulos (b, h, y_centro): área total A, altura del centroide ȳ, momento de inercia I
    respecto al eje centroidal horizontal (Σ b h³/12 + A d²), la suma sin Steiner y el detalle por rectángulo."""
    rects = list(rectangulos)
    if not rects:
        raise ValueError("escribe al menos un rectángulo.")
    datos = []
    for k, rect in enumerate(rects, 1):
        if len(rect) != 3:
            raise ValueError(f"el rectángulo {k} necesita tres números: b, h, y_centro.")
        b, h, yc = (exacto(v, f"el dato del rectángulo {k}") for v in rect)
        if b <= 0 or h <= 0:
            raise ValueError(f"en el rectángulo {k} la base y la altura deben ser positivas.")
        datos.append((b, h, yc))
    A = sum(b * h for b, h, _ in datos)
    yb = sum(b * h * yc for b, h, yc in datos) / A
    detalle = [(b * h, yc - yb, b * h**3 / 12, b * h * (yc - yb)**2) for b, h, yc in datos]
    I_propias = sum(d[2] for d in detalle)
    I = I_propias + sum(d[3] for d in detalle)
    return {"A": A, "y_bar": yb, "I": I, "I_sin_steiner": I_propias, "detalle": detalle}


def _enciman(rects):
    """True si dos rectángulos ocupan la misma franja de alturas (centrados en el mismo eje, se enciman)."""
    franjas = sorted((yc - h / 2, yc + h / 2) for _, h, yc in rects)
    return any(s2 < e1 - 1e-12 for (_, e1), (s2, _) in zip(franjas, franjas[1:]))


def calculadora_seccion(r1_txt, r2_txt, r3_txt):
    try:
        rects = [lista_numeros(t, f"el rectángulo {k}") for k, t in enumerate((r1_txt, r2_txt, r3_txt), 1) if t.strip()]
        r = seccion_compuesta(rects)
    except ERRORES as err:
        print("Revisa la entrada:", explica(err)); return
    print(f"({CIFRAS} cifras significativas, redondeadas; longitudes en cm)")
    print(" rect.   A_i (cm²)    d_i (cm)   b h³/12 (cm⁴)   A_i d_i² (cm⁴)")
    for k, (Ai, di, Ii, Ad2) in enumerate(r["detalle"], 1):
        print(f"  {k}   {cifras(Ai, CIFRAS):>10}  {cifras(di, CIFRAS):>10}  {cifras(Ii, CIFRAS):>14}  {cifras(Ad2, CIFRAS):>14}")
    print(f"A = {cifras(r['A'], CIFRAS)} cm²     ȳ = {cifras(r['y_bar'], CIFRAS)} cm (desde y = 0)")
    print(f"I respecto al eje centroidal = {cifras(r['I'], CIFRAS)} cm⁴")
    if r["I"] - r["I_sin_steiner"] > 0:
        falta = 100 * (1 - r["I_sin_steiner"] / r["I"])
        print(f"sin el término A d² saldría {cifras(r['I_sin_steiner'], CIFRAS)} cm⁴, un {cifras(falta, 3)} % menos "
              "(porcentaje con 3 cifras): ese resultado es incorrecto")
    if _enciman(rects):
        print("Aviso: dos rectángulos ocupan la misma franja de alturas. Si en la pieza se tocan, el área común "
              "se está contando dos veces.")
    fig, ax = plt.subplots(figsize=(4.5, 3.6))
    for (b, h, yc), achurado in zip(rects, ("///", "\\\\\\", "xxx")):
        ax.add_patch(plt.Rectangle((-b / 2, yc - h / 2), b, h, facecolor="white", edgecolor="black", hatch=achurado))
    yb = float(r["y_bar"])
    ancho = max(b for b, _, _ in rects)
    ax.axhline(yb, color="black", ls="--", lw=1)
    rotula(ax, ancho / 2, yb, f"eje centroidal, ȳ = {cifras(yb, CIFRAS)} cm", xytext=(6, 8))
    ax.autoscale(); ax.set_aspect("equal"); ax.margins(0.15)
    ax.set_xlabel("cm"); ax.set_ylabel("y (cm)"); plt.show()


_ancho = {"description_width": "initial"}
widgets.interact(calculadora_seccion,
    r1_txt=widgets.Text(value="12, 2, 11", description="rect. 1 (b, h, y):", continuous_update=False, style=_ancho),
    r2_txt=widgets.Text(value="2, 10, 5", description="rect. 2 (b, h, y):", continuous_update=False, style=_ancho),
    r3_txt=widgets.Text(value="", description="rect. 3 (b, h, y):", continuous_update=False, style=_ancho));

# %%
# Casos de prueba de la sección 6.4 (resultado conocido)
def _seccion(rects):
    return {k: float(v) for k, v in seccion_compuesta(rects).items() if k != "detalle"}


PRUEBAS_4 = [
    ("región bajo 4 - x² en [0, 2]: A = 16/3, x̄ = 3/4, ȳ = 8/5",
     lambda: (lambda r: (r["A"], r["x_bar"], r["y_bar"]) == (sp.Rational(16, 3), sp.Rational(3, 4), sp.Rational(8, 5)))(
         region(4 - x**2, 0, 0, 2))),
    ("rectángulo 4 × 6 (f = 6, g = 0 en [0, 4]): I_x = 288 e I centroidal = 72 = b h³/12",
     lambda: (lambda r: r["I_x"] == 288 and r["I_c"] == 72)(region(6 + 0 * x, 0, 0, 4))),
    ("triángulo de base 4 y altura 6 (bajo y = 3x/2 en [0, 4]): centroide (8/3, 2)",
     lambda: (lambda r: (r["x_bar"], r["y_bar"]) == (sp.Rational(8, 3), 2))(region(3 * x / 2, 0, 0, 4))),
    ("región con g arriba de f se rechaza", lambda: _rechaza(lambda: region(0, 4 - x**2, 0, 2))),
    ("sección T (patín 12 × 2 en y = 11, alma 2 × 10 en y = 5): ȳ = 8.2727, I = 567.39 cm⁴, sin Steiner 174.67",
     lambda: (lambda r: round(r["y_bar"], 4) == 8.2727 and round(r["I"], 2) == 567.39
              and round(r["I_sin_steiner"], 2) == 174.67)(_seccion([(12, 2, 11), (2, 10, 5)]))),
    ("sección I (patines 19 × 1 en y = ±14.5, alma 0.6 × 28): I = 9090.27 cm⁴, A = 54.8 cm²",
     lambda: (lambda r: round(r["I"], 2) == 9090.27 and cerca(r["A"], 54.8))(
         _seccion([(19, 1, 14.5), (19, 1, -14.5), (0.6, 28, 0)]))),
    ("rectángulo 4 × 30: I = 9000 cm⁴", lambda: seccion_compuesta([(4, 30, 15)])["I"] == 9000),
    ("las dos calculadoras coinciden: rectángulo 4 × 6 da 72 en ambas",
     lambda: seccion_compuesta([(4, 6, 3)])["I"] == region(6 + 0 * x, 0, 0, 4)["I_c"] == 72),
    ("base negativa se rechaza", lambda: _rechaza(lambda: seccion_compuesta([(-2, 10, 5)]))),
]
for nombre, prueba in PRUEBAS_4:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 6.4).** Un triángulo rectángulo tiene base 6 sobre el eje $x$ y altura 3. Calcula a mano la altura $\bar y$ de su centroide y escribe el número.

# %%
mi_valor = None     # escribe un número

if mi_valor is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    try:
        mio = a_numero(mi_valor)
    except ValueError as err:
        print("Revisa tu número:", err)
    else:
        ref = region(x / 2, 0, 0, 6)["y_bar"]          # hipotenusa y = x/2 de (0, 0) a (6, 3)
        print("La calculadora da: ȳ =", muestra(ref))
        print("coinciden a 3 cifras" if coincide(mio, ref) else
              "NO coinciden: con f = x/2 y g = 0, ȳ = (1/A) ∫ f²/2 dx; el centroide de un triángulo queda a h/3 de la base")
