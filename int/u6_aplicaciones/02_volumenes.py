# ID: INT-U6-NB02
# Notebook: int/u6_aplicaciones.ipynb · sección 6.2 volúmenes de sólidos de revolución
# Repositorio: int/u6_aplicaciones/02_volumenes.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6.2 Volúmenes de sólidos de revolución
#
# $$V_{\text{arandelas}}=\pi\int_a^b\Big(f(x)^2-g(x)^2\Big)\,dx\quad(\text{giro alrededor del eje }x),\qquad V_{\text{cascarones}}=2\pi\int_a^b x\,\big(f(x)-g(x)\big)\,dx\quad(\text{giro alrededor del eje }y).$$
#
# Con arandelas, el radio exterior es la curva más alejada del eje y el interior la más cercana; si $g=0$ quedan discos. Con cascarones, $x$ es el radio de cada cascarón y $f-g$ su altura, así que la región debe quedar a la derecha del eje $y$ ($a\ge0$).
#
# En el ejemplo exacto, el cono de radio $r$ y altura $h$ sale de girar $y=\frac{r}{h}x$ en $[0,h]$ alrededor del eje $x$.

# %%
r_c, h_c = sp.symbols("r h", positive=True)
cono = sp.pi * sp.integrate((r_c / h_c * x)**2, (x, 0, h_c))
print("V del cono =", sp.factor(cono), "   con r = 3, h = 4:", cono.subs({r_c: 3, h_c: 4}))

# %% [markdown]
# **Calculadora de volúmenes.** Escribe $f$ y $g$ en cualquier orden (para discos, $g=0$), los extremos y el método. Con arandelas, en cada $x$ el radio exterior es el mayor de $|f|$ y $|g|$ y el interior el menor; la calculadora rechaza una región que quede a los dos lados del eje de giro, porque al girar sus dos partes se enciman. La gráfica dibuja la región en gris y su reflejo respecto al eje de giro, que es la silueta del sólido.

# %%
def _suma_de_magnitudes(h, a, b, integ):
    """Σ |∫ h| sobre las partes de [a, b] donde h no cambia de signo (∫ |h| dx sin integrar un valor absoluto)."""
    puntos, _ = cruces(h, a, b)
    integ.solo_numerica = integ.solo_numerica or len(puntos) > 20
    return simplifica(sum(abs(integ(h, p, q)) for p, q in zip([a, *puntos], [*puntos, b])))


def volumen(f, g, a, b, metodo="arandelas"):
    """Volumen del sólido que genera la región entre y = f(x) y y = g(x) en [a, b] (f y g en cualquier orden).
    metodo "arandelas": giro alrededor del eje x; en cada x el radio exterior es el mayor de |f| y |g| y el interior
    el menor (discos si g = 0); se rechaza una región que quede a los dos lados del eje x.
    metodo "cascarones": giro alrededor del eje y; radio x y altura |f - g|; exige a >= 0.
    Devuelve un dict con V, el integrando (sin el valor absoluto) y el método de integración."""
    if metodo not in ("arandelas", "cascarones"):
        raise ValueError('el método es "arandelas" o "cascarones".')
    f, g = sp.sympify(f), sp.sympify(g)
    a, b = exacto(a, "a"), exacto(b, "b")
    if not a < b:
        raise ValueError("a debe ser menor que b.")
    revisa_intervalo(f, a, b, "f"); revisa_intervalo(g, a, b, "g")
    if simplifica(f - g) == 0:
        raise ValueError("f y g son la misma función: no encierran ninguna región (el volumen vale 0).")
    integ = Integrador()
    if metodo == "arandelas":
        xs = malla_interior(a, b)
        F, G = numerica(f)(xs), numerica(g)(xs)
        tol = 1e-9 * max(1.0, float(np.max(np.abs(F))), float(np.max(np.abs(G))))
        lados = ((F > tol) & (G < -tol)) | ((F < -tol) & (G > tol))
        if np.any(lados):
            raise ValueError(f"en x ≈ {cifras(xs[int(np.argmax(lados))], 4)} la región queda a los dos lados del eje x "
                             "(f y g tienen signos opuestos): al girar, la parte de arriba y la de abajo se enciman. "
                             "Calcula cada parte por separado; el sólido lo genera la más alejada del eje.")
        integrando = sp.pi * (f**2 - g**2)
    else:
        if a < 0:
            raise ValueError("con cascarones el radio es x y no puede ser negativo: usa a ≥ 0, con la región a la "
                             "derecha del eje y.")
        integrando = 2 * sp.pi * x * (f - g)
    V = _suma_de_magnitudes(sp.expand(integrando), a, b, integ)
    return {"V": V, "integrando": integrando, "metodo": integ.metodo}


def calculadora_volumen(f_txt, g_txt, a_txt, b_txt, metodo):
    try:
        f, g = parsear(f_txt), parsear(g_txt)
        a, b = exacto(a_txt, "a"), exacto(b_txt, "b")
        r = volumen(f, g, a, b, metodo)
    except ERRORES as err:
        print("Revisa la entrada:", explica(err)); return
    eje = "eje x" if metodo == "arandelas" else "eje y"
    print(f"({CIFRAS} cifras significativas, redondeadas; si x y y están en metros, el volumen sale en m³)")
    print(f"giro alrededor del {eje}:  V = ∫ de {punto(a)} a {punto(b)} de  |{r['integrando']}|  dx")
    print(f"V = {muestra(r['V'])} u³   ({COMO[r['metodo']]})")
    try:
        xs = np.linspace(float(a), float(b), 400)
        yf, yg = para_graficar(f, xs), para_graficar(g, xs)
        fig, ax = plt.subplots(figsize=(5.5, 3.2))
        k1, k2 = len(xs) * 2 // 5, len(xs) * 7 // 10
        if metodo == "arandelas":
            for s_ in (1, -1):
                ax.fill_between(xs, s_ * yg, s_ * yf, color="0.8" if s_ == 1 else "0.92")
                ax.plot(xs, s_ * yf, color="black", ls="-"); ax.plot(xs, s_ * yg, color="black", ls="--")
            ax.axhline(0, color="black", lw=0.8, ls="-.")
            rotula(ax, xs[-1], 0, "eje de giro", ha="right", xytext=(0, 7))
        else:
            for s_ in (1, -1):
                ax.fill_between(s_ * xs, yg, yf, color="0.8" if s_ == 1 else "0.92")
                ax.plot(s_ * xs, yf, color="black", ls="-"); ax.plot(s_ * xs, yg, color="black", ls="--")
            ax.axvline(0, color="black", lw=0.8, ls="-.")
            rotula(ax, 0, float(np.nanmax(np.concatenate((yf, yg)))), "eje de giro", xytext=(4, 6))
        rotula(ax, xs[k1], yf[k1], "f", xytext=(0, 9), ha="center")
        rotula(ax, xs[k2], yg[k2], "g", xytext=(0, -9), ha="center")
        ax.set_xlabel("x"); ax.set_title(f"región (gris) y su reflejo: silueta del sólido, {metodo}", fontsize=9)
        plt.show()
    except ERRORES:
        print("(no pude dibujar la gráfica con estos datos)")


widgets.interact(calculadora_volumen,
    f_txt=widgets.Text(value="sqrt(x)", description="f(x) =", continuous_update=False),
    g_txt=widgets.Text(value="x/2", description="g(x) =", continuous_update=False),
    a_txt=widgets.Text(value="0", description="a =", continuous_update=False),
    b_txt=widgets.Text(value="4", description="b =", continuous_update=False),
    metodo=widgets.Dropdown(options=[("arandelas o discos (eje x)", "arandelas"), ("cascarones (eje y)", "cascarones")],
                            value="arandelas", description="método"));

# %%
# Casos de prueba de la sección 6.2 (resultado conocido)
PRUEBAS_2 = [
    ("discos de y = 3x/4 en [0, 4]: cono de 12π", lambda: volumen(3 * x / 4, 0, 0, 4, "arandelas")["V"] == 12 * sp.pi),
    ("arandelas con √x y x/2 en [0, 4], en cualquier orden: 8π/3",
     lambda: volumen(sp.sqrt(x), x / 2, 0, 4, "arandelas")["V"] == volumen(x / 2, sp.sqrt(x), 0, 4, "arandelas")["V"]
     == sp.Rational(8, 3) * sp.pi),
    ("cascarones con 4x - x² en [0, 4]: 128π/3",
     lambda: volumen(4 * x - x**2, 0, 0, 4, "cascarones")["V"] == sp.Rational(128, 3) * sp.pi),
    ("región bajo el eje x entre -x y -2x en [0, 1]: radios x y 2x, V = π",
     lambda: volumen(-x, -2 * x, 0, 1, "arandelas")["V"] == sp.pi),
    ("discos de y = x en [-1, 1]: dos conos, V = 2π/3", lambda: volumen(x, 0, -1, 1, "arandelas")["V"] == 2 * sp.pi / 3),
    ("se rechazan una región a los dos lados del eje x (1 y -1) y cascarones con a < 0",
     lambda: _rechaza(lambda: volumen(1, -1, 0, 1, "arandelas"))
     and _rechaza(lambda: volumen(1 - x**2, 0, -1, 1, "cascarones"))),
]
for nombre, prueba in PRUEBAS_2:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 6.2).** Calcula a mano, por discos, el volumen de la esfera de radio 1 que se forma al girar $y=\sqrt{1-x^2}$ en $[-1,1]$ alrededor del eje $x$, y escribe el número.

# %%
mi_valor = None     # escribe un número (también sirve "4*pi/3")

if mi_valor is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    try:
        mio = a_numero(mi_valor)
    except ValueError as err:
        print("Revisa tu número:", err)
    else:
        ref = volumen(sp.sqrt(1 - x**2), 0, -1, 1, "arandelas")["V"]
        print("La calculadora da:", muestra(ref), "u³")
        print("coinciden a 3 cifras" if coincide(mio, ref) else
              "NO coinciden: el radio del disco es √(1 - x²), así que su área es π(1 - x²); integra de -1 a 1")
