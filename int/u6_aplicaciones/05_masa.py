# ID: INT-U6-NB05
# Notebook: int/u6_aplicaciones.ipynb · sección 6.5 momentos de masa
# Repositorio: int/u6_aplicaciones/05_masa.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6.5 Momentos de masa
#
# Para una barra sobre $[0,L]$ con densidad lineal $\lambda(x)$ en kg/m:
#
# $$m=\int_0^L\lambda(x)\,dx,\qquad \bar x=\frac1m\int_0^Lx\,\lambda(x)\,dx,\qquad I_O=\int_0^Lx^2\,\lambda(x)\,dx,\qquad I_G=I_O-m\,\bar x^{\,2}.$$
#
# $I_O$ es el momento de inercia respecto a un eje perpendicular a la barra que pasa por el extremo $x=0$, e $I_G$ el que pasa por el centro de masa. Una masa puntual $m_p$ en la punta $x=L$ suma $m_p$ a la masa, $m_pL$ al primer momento y $m_pL^2$ a $I_O$.
#
# El ejemplo exacto es una barra de $L=\tfrac12$ m con $\lambda(x)=4-4x$ kg/m, que se adelgaza hacia la punta.

# %%
lam_ej, L_ej5 = 4 - 4 * x, sp.Rational(1, 2)
m_ej = sp.integrate(lam_ej, (x, 0, L_ej5))
print("m =", m_ej, "kg   x̄ =", sp.integrate(x * lam_ej, (x, 0, L_ej5)) / m_ej, "m   I_O =",
      sp.integrate(x**2 * lam_ej, (x, 0, L_ej5)), "kg·m²")

# %% [markdown]
# **Calculadora de una barra.** Escribe $\lambda(x)$ en kg/m, la longitud $L$ en m y, si la hay, la masa puntual en la punta en kg. La calculadora da la masa, el centro de masa y los momentos de inercia respecto al extremo y al centro de masa. En la gráfica, el área bajo $\lambda$ es la masa de la barra.

# %%
def barra(lam, L, m_punta=0):
    """Barra sobre [0, L] con densidad lineal lam(x) en kg/m y una masa puntual opcional m_punta en x = L.
    Devuelve un dict con m (total), x_bar, I_O (respecto a x = 0), I_G (respecto al centro de masa) y m_barra."""
    lam = sp.sympify(lam)
    L, mp = exacto(L, "L"), exacto(m_punta, "la masa en la punta")
    if L <= 0:
        raise ValueError("la longitud L debe ser positiva.")
    if mp < 0:
        raise ValueError("la masa en la punta no puede ser negativa (escribe 0 si no hay).")
    revisa_intervalo(lam, 0, L, "λ")
    if lam.has(x):
        valores = numerica(lam)(malla_interior(0, L))
        if np.any(valores < -1e-12):
            raise ValueError(f"λ es negativa cerca de x ≈ {cifras(malla_interior(0, L)[int(np.argmin(valores))], 4)} m: "
                             "una densidad no puede ser negativa. Revisa la fórmula o la longitud.")
    elif lam < 0:
        raise ValueError("una densidad no puede ser negativa.")
    integ = Integrador()
    def _int(h):
        return integ(sp.expand(h), 0, L)
    m_b = _int(lam)
    m = m_b + mp
    if float(m) <= 0:
        raise ValueError("la masa total es cero: escribe una densidad positiva o una masa en la punta.")
    xb = simplifica((_int(x * lam) + mp * L) / m)
    IO = simplifica(_int(x**2 * lam) + mp * L**2)
    IG = simplifica(IO - m * xb**2)
    return {"m": simplifica(m), "x_bar": xb, "I_O": IO, "I_G": IG, "m_barra": m_b,
            "metodo": integ.metodo}


def calculadora_barra(lam_txt, L, m_punta):
    try:
        lam = parsear(lam_txt)
        r = barra(lam, L, m_punta)
    except ERRORES as err:
        print("Revisa la entrada:", explica(err)); return
    print(f"({CIFRAS} cifras significativas, redondeadas; {COMO[r['metodo']]})")
    print(f"masa de la barra = {muestra(r['m_barra'])} kg" + (f"   masa total = {muestra(r['m'])} kg" if m_punta else ""))
    print(f"x̄   = {muestra(r['x_bar'])} m   (centro de masa, medido desde x = 0)")
    print(f"I_O = {muestra(r['I_O'])} kg·m²   (eje en x = 0)")
    print(f"I_G = {muestra(r['I_G'])} kg·m²   (eje en el centro de masa, I_O - m x̄²)")
    try:
        xs = np.linspace(0, float(L), 300)
        ys = para_graficar(lam, xs)
        fig, ax = plt.subplots(figsize=(5.5, 3))
        ax.fill_between(xs, 0, ys, color="0.85"); ax.plot(xs, ys, color="black")
        rotula(ax, xs[len(xs) // 3], ys[len(ys) // 3], "λ(x)", xytext=(4, 8))
        xb = float(r["x_bar"])
        ax.axvline(xb, color="black", ls="--", lw=1)
        rotula(ax, xb, float(np.nanmax(ys)) * 0.5, f"x̄ = {cifras(xb, CIFRAS)} m")
        if m_punta:
            ax.plot([float(L)], [0], marker="s", color="black", ms=9, clip_on=False)
            rotula(ax, float(L), 0, f"{m_punta:g} kg", xytext=(-6, 12), ha="right")
        ax.set_xlabel("x (m)"); ax.set_ylabel("λ (kg/m)"); ax.set_ylim(bottom=0)
        ax.set_title("el área gris es la masa de la barra", fontsize=9); plt.show()
    except ERRORES:
        print("(no pude dibujar la gráfica con estos datos)")


widgets.interact(calculadora_barra,
    lam_txt=widgets.Text(value="4 - 4x", description="λ(x) =", continuous_update=False),
    L=widgets.FloatText(value=0.5, description="L (m)"),
    m_punta=widgets.FloatText(value=0.5, description="punta (kg)"));

# %%
# Casos de prueba de la sección 6.5 (resultado conocido)
PRUEBAS_5 = [
    ("λ = 2, L = 3: m = 6, x̄ = 1.5, I_O = 18",
     lambda: (lambda r: (r["m"], r["x_bar"], r["I_O"]) == (6, sp.Rational(3, 2), 18))(barra(2, 3))),
    ("λ = 4 - 4x, L = 1/2: m = 3/2, x̄ = 2/9, I_O = 5/48",
     lambda: (lambda r: (r["m"], r["x_bar"], r["I_O"]) == (sp.Rational(3, 2), sp.Rational(2, 9), sp.Rational(5, 48)))(
         barra(4 - 4 * x, sp.Rational(1, 2)))),
    ("λ = x, L = 1: m = 1/2, x̄ = 2/3, I_O = 1/4",
     lambda: (lambda r: (r["m"], r["x_bar"], r["I_O"]) == (sp.Rational(1, 2), sp.Rational(2, 3), sp.Rational(1, 4)))(
         barra(x, 1))),
    ("λ = 4 - 4x, L = 1/2 con 0.5 kg en la punta: I_O = 0.229167",
     lambda: cifras(barra(4 - 4 * x, 0.5, 0.5)["I_O"], 6) == "0.229167"),
    ("barra uniforme de 6 kg y 3 m: I_G = m L²/12 = 4.5", lambda: barra(2, 3)["I_G"] == sp.Rational(9, 2)),
    ("densidad negativa se rechaza", lambda: _rechaza(lambda: barra(1 - 2 * x, 1))),
]
for nombre, prueba in PRUEBAS_5:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 6.5).** Una barra uniforme de 2 kg mide 1 m. Calcula a mano su momento de inercia $I_O$ respecto a un eje perpendicular que pasa por un extremo y escribe el número en kg·m².

# %%
mi_valor = None     # escribe un número (también sirve "2/3")

if mi_valor is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    try:
        mio = a_numero(mi_valor)
    except ValueError as err:
        print("Revisa tu número:", err)
    else:
        ref = barra(2, 1)["I_O"]          # λ = 2 kg/m en [0, 1]
        print("La calculadora da: I_O =", muestra(ref), "kg·m²")
        print("coinciden a 3 cifras" if coincide(mio, ref) else
              "NO coinciden: con λ = m/L = 2 kg/m, I_O = ∫_0^1 2x² dx; recuerda que el resultado general es m L²/3")
