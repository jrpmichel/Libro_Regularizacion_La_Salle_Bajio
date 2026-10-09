# ID: INT-U6-NB09
# Notebook: int/u6_aplicaciones.ipynb · sección 6.9 carga, energía y valor eficaz de una señal
# Repositorio: int/u6_aplicaciones/09_carga.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6.9 Carga, energía y valor eficaz de una señal
#
# $$q=\int_{t_0}^{t_1}i(t)\,dt\ \ (\text{C}=\text{A·s}),\qquad E=\int_{t_0}^{t_1}v(t)\,i(t)\,dt\ \ (\text{J}),\qquad f_{\text{ef}}=\sqrt{\frac1T\int_0^Tf(t)^2\,dt}.$$
#
# El valor eficaz (RMS) es el valor constante que disiparía la misma potencia en una resistencia. Para un sensor que despierta, mide, transmite y vuelve a dormir cada periodo $T_p$, la corriente media y la vida de la batería son
#
# $$\bar I=\frac{q_{\text{activo}}+i_{\text{reposo}}\,(T_p-t_{\text{activo}})}{T_p},\qquad \text{vida}=\frac{\text{capacidad}\cdot\text{fracción utilizable}}{\bar I},\qquad 1\ \text{mAh}=3.6\ \text{C}.$$
#
# El ejemplo exacto integra la corriente de arranque de un sensor, $i(t)=3+9e^{-t/0.4}$ mA, durante 1.5 s.

# %%
t_c = sp.symbols("t", nonnegative=True)
q_ej = sp.integrate(3 + 9 * sp.exp(-t_c / sp.Rational(2, 5)), (t_c, 0, sp.Rational(3, 2)))
print("q =", q_ej, "mC ≈", cifras(q_ej, 7), "mC")

# %% [markdown]
# **Calculadora de carga y energía.** Escribe $i(t)$ y elige si está en A o en mA; escribe también el voltaje $v(t)$ en V (puede ser constante, o vacío si no lo necesitas) y el intervalo en s (los extremos aceptan fracciones como `1/60`). La calculadora da la carga, la corriente media, la corriente eficaz y la energía entregada.

# %%
def carga(i, t0, t1, con_metodo=False):
    """q = ∫_{t0}^{t1} i(t) dt: en C si i está en A (en mC si está en mA). Con con_metodo=True devuelve (q, método)."""
    i = sp.sympify(i)
    t0, t1 = exacto(t0, "t0"), exacto(t1, "t1")
    if not t0 < t1:
        raise ValueError("t0 debe ser menor que t1.")
    revisa_intervalo(i, t0, t1, "i", var="t")
    q, como = integra(i, t0, t1)
    return (q, como) if con_metodo else q


def energia(v, i, t0, t1):
    """E = ∫_{t0}^{t1} v(t) i(t) dt: en J si v está en V e i en A."""
    return carga(sp.sympify(v) * sp.sympify(i), t0, t1)


def valor_eficaz(f, T, t0=0):
    """Valor eficaz (RMS) de f en un periodo [t0, t0 + T]: √((1/T) ∫ f² dt)."""
    f = sp.sympify(f)
    T, t0 = exacto(T, "T"), exacto(t0, "t0")
    if T <= 0:
        raise ValueError("el periodo T debe ser positivo.")
    return simplifica(sp.sqrt(carga(f**2, t0, t0 + T) / T))


def vida_bateria(q_activo, t_activo, i_reposo, periodo, capacidad_mAh, margen):
    """Corriente media (A) y vida de la batería (años de 365 días) para un ciclo que gasta q_activo (C) durante
    t_activo (s) y luego i_reposo (A) hasta completar el periodo (s). margen es la fracción utilizable (0.8 = 80 %)."""
    q_activo, t_activo, i_reposo = float(q_activo), float(t_activo), float(i_reposo)
    periodo, capacidad_mAh, margen = float(periodo), float(capacidad_mAh), float(margen)
    if q_activo < 0 or i_reposo < 0:
        raise ValueError("la carga activa y la corriente de reposo no pueden ser negativas.")
    if periodo <= 0 or not 0 <= t_activo <= periodo:
        raise ValueError("el periodo debe ser positivo y el tiempo activo debe caber dentro de él.")
    if capacidad_mAh <= 0:
        raise ValueError("la capacidad debe ser positiva.")
    if not 0 < margen <= 1:
        raise ValueError("el margen es la fracción utilizable de la capacidad: un número entre 0 y 1 (0.8 = 80 %).")
    I_media = (q_activo + i_reposo * (periodo - t_activo)) / periodo
    if I_media <= 0:
        raise ValueError("la corriente media es cero: la batería no se gastaría.")
    vida_s = capacidad_mAh * 3.6 * margen / I_media
    return I_media, vida_s / (365 * 24 * 3600)


def calculadora_carga(i_txt, unidad, v_txt, t0_txt, t1_txt):
    try:
        i = parsear(i_txt)
        t0, t1 = exacto(t0_txt, "t0"), exacto(t1_txt, "t1")
        q, como = carga(i, t0, t1, con_metodo=True)
        T = float(t1 - t0)
        ts = np.linspace(float(t0), float(t1), 500)
        ys = para_graficar(i, ts)
    except ERRORES as err:
        print("Revisa la entrada:", explica(err)); return
    uq, ue = ("C", "J") if unidad == "A" else ("mC", "mJ")
    finitos = ys[np.isfinite(ys)]
    escala = float(np.max(np.abs(finitos))) * T if finitos.size else abs(float(q))
    nula = abs(float(q)) <= 1e-12 * escala            # cancelación numérica: q es 0
    print(f"({CIFRAS} cifras significativas, redondeadas; t en s; {COMO[como]})")
    print(f"carga q = {'0' if nula else muestra(q)} {uq}")
    print(f"corriente media = {'0' if nula else cifras(float(q) / T, CIFRAS)} {unidad}")
    try:
        print(f"corriente eficaz = {cifras(valor_eficaz(i, t1 - t0, t0), CIFRAS)} {unidad}")
    except ERRORES as err:
        print("corriente eficaz: no existe, porque la integral de i² diverge" if "diverge" in str(err)
              else f"corriente eficaz: {explica(err)}")
    if v_txt.strip():
        try:
            print(f"energía E = {muestra(energia(parsear(v_txt), i, t0, t1))} {ue}")
        except ERRORES as err:
            print("energía:", explica(err))
    try:
        fig, ax = plt.subplots(figsize=(5.5, 3))
        ax.fill_between(ts, 0, ys, color="0.85"); ax.plot(ts, ys, color="black")
        im = 0.0 if nula else float(q) / T
        ax.axhline(im, color="black", ls="--", lw=1); rotula(ax, ts[-1], im, "media", xytext=(-30, 7))
        ax.set_xlabel("t (s)"); ax.set_ylabel(f"i ({unidad})"); ax.set_ylim(bottom=min(0, float(np.nanmin(ys))))
        ax.set_title(f"el área gris es la carga q, en {uq}", fontsize=9); plt.show()
    except ERRORES:
        print("(no pude dibujar la gráfica con estos datos)")


widgets.interact(calculadora_carga,
    i_txt=widgets.Text(value="3 + 9*exp(-t/0.4)", description="i(t) =", continuous_update=False),
    unidad=widgets.Dropdown(options=["A", "mA"], value="mA", description="i en"),
    v_txt=widgets.Text(value="3.3", description="v(t) (V) =", continuous_update=False),
    t0_txt=widgets.Text(value="0", description="t0 (s)", continuous_update=False),
    t1_txt=widgets.Text(value="1.5", description="t1 (s)", continuous_update=False));

# %% [markdown]
# **Calculadora de la vida de una batería.** Escribe la carga que gasta el sensor en cada ciclo activo (en mC; por ejemplo, los 8.015 mC del arranque más 40 mA durante 0.12 s de transmisión dan 12.815 mC en 1.62 s), la corriente de reposo, el periodo y la batería. La gráfica muestra cómo cambia la vida con el periodo.

# %%
def calculadora_bateria(q_mC, t_activo, i_reposo_uA, periodo, capacidad_mAh, util_pct):
    try:
        I, anios = vida_bateria(q_mC / 1000, t_activo, i_reposo_uA * 1e-6, periodo, capacidad_mAh, util_pct / 100)
    except ERRORES as err:
        print("Revisa la entrada:", explica(err)); return
    print(f"({CIFRAS} cifras significativas, redondeadas; años de 365 días)")
    print(f"corriente media = {cifras(I * 1e6, CIFRAS)} µA")
    print(f"vida = {cifras(anios, CIFRAS)} años = {cifras(anios * 365, CIFRAS)} días")
    try:
        periodos = np.linspace(max(t_activo * 1.05, 0.1 * periodo), 5 * periodo, 200)
        vidas = [vida_bateria(q_mC / 1000, t_activo, i_reposo_uA * 1e-6, p, capacidad_mAh, util_pct / 100)[1]
                 for p in periodos]
        fig, ax = plt.subplots(figsize=(5.5, 3))
        ax.plot(periodos, vidas, color="black"); ax.plot([periodo], [anios], marker="o", color="black")
        rotula(ax, periodo, anios, f"tu diseño: {cifras(anios, CIFRAS)} años", xytext=(6, -10))
        ax.set_xlabel("periodo entre mediciones (s)"); ax.set_ylabel("vida (años)"); ax.set_ylim(bottom=0)
        ax.set_title("vida de la batería contra el periodo", fontsize=9); plt.show()
    except ERRORES:
        print("(no pude dibujar la gráfica con estos datos)")


_ancho = {"description_width": "initial"}
widgets.interact(calculadora_bateria,
    q_mC=widgets.FloatText(value=12.815336, description="carga activa (mC)", style=_ancho),
    t_activo=widgets.FloatText(value=1.62, description="tiempo activo (s)", style=_ancho),
    i_reposo_uA=widgets.FloatText(value=6, description="corriente de reposo (µA)", style=_ancho),
    periodo=widgets.FloatText(value=120, description="periodo (s)", style=_ancho),
    capacidad_mAh=widgets.FloatText(value=2600, description="capacidad (mAh)", style=_ancho),
    util_pct=widgets.IntSlider(value=80, min=10, max=100, description="utilizable (%)", style=_ancho,
                               continuous_update=False));

# %%
# Casos de prueba de la sección 6.9 (resultado conocido)
_pwm = sp.Piecewise((1, x < sp.Rational(1, 4)), (0, True))

PRUEBAS_9 = [
    ("i = 2 A de 0 a 5 s: 10 C", lambda: carga(2, 0, 5) == 10),
    ("i = 3 + 9e^(-t/0.4) mA de 0 a 1.5 s: 8.015336 mC",
     lambda: cifras(carga(3 + 9 * sp.exp(-x / sp.Rational(2, 5)), 0, 1.5), 7) == "8.015336"),
    ("i = 0.5e^(-50t) A de 0 a 0.05 s: 0.0091792 C",
     lambda: cifras(carga(sp.exp(-50 * x) / 2, 0, 0.05), 5) == "0.0091792"),
    ("v = 3.3 V e i = 2 A durante 5 s: 33 J", lambda: energia(sp.Rational("3.3"), 2, 0, 5) == 33),
    ("valor eficaz de la sierra v = t en [0, 1): 1/√3",
     lambda: sp.simplify(valor_eficaz(x, 1) - 1 / sp.sqrt(3)) == 0),
    ("valor eficaz de un PWM de amplitud 1 y ciclo de trabajo 0.25: 0.5", lambda: valor_eficaz(_pwm, 1) == sp.Rational(1, 2)),
    ("el mismo PWM escrito con escalón: 0.5", lambda: valor_eficaz(parsear("escalon(0.25 - t)"), 1) == sp.Rational(1, 2)),
    ("vida del sensor (12.815336 mC en 1.62 s, 6 µA, 120 s, 2600 mAh, 80 %): 112.71 µA y 2.1066 años",
     lambda: (lambda r: cifras(r[0] * 1e6, 5) == "112.71" and cifras(r[1], 5) == "2.1066")(
         vida_bateria(12.815336e-3, 1.62, 6e-6, 120, 2600, 0.8))),
    ("3.6 C cada hora, sin reposo (1 mA de media) con 2600 mAh al 100 %: 2600 h",
     lambda: (lambda r: cerca(r[0], 1e-3, 1e-15) and cerca(r[1], 2600 / 8760, 1e-12))(
         vida_bateria(3.6, 3600, 0, 3600, 2600, 1))),
    ("1 mA constante con 1000 mAh: 1000 h = 1000/8760 años",
     lambda: cerca(vida_bateria(0, 0, 1e-3, 1, 1000, 1)[1], 1000 / 8760, 1e-12)),
    ("margen mayor que 1 se rechaza", lambda: _rechaza(lambda: vida_bateria(0.01, 1, 6e-6, 120, 2600, 80))),
    ("tiempo activo mayor que el periodo se rechaza", lambda: _rechaza(lambda: vida_bateria(0.01, 200, 6e-6, 120, 2600, 0.8))),
]
for nombre, prueba in PRUEBAS_9:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 6.9).** Un capacitor se descarga con $i(t)=0.5\,e^{-50t}$ A. Calcula a mano la carga que pasa entre $t=0$ y $t=0.1$ s y escribe el número en C.

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
        ref = carga(sp.exp(-50 * x) / 2, 0, sp.Rational(1, 10))
        print("La calculadora da: q =", muestra(ref), "C")
        print("coinciden a 3 cifras" if coincide(mio, ref) else
              "NO coinciden: la antiderivada de 0.5 e^(-50t) es -0.01 e^(-50t); evalúa entre 0 y 0.1")
