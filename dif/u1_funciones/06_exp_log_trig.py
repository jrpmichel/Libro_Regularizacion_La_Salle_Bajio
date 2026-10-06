# ID: DIF-U1-NB06
# Notebook: dif/u1_funciones.ipynb · sección 1.6 exponencial, logaritmo y trigonométricas
# Repositorio: dif/u1_funciones/06_exp_log_trig.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6. Funciones exponencial, logarítmica y trigonométricas
#
# Tres calculadoras. La primera da seno, coseno y tangente con la **unidad del ángulo visible** (grados o radianes). La segunda resuelve $a^x=c$ con logaritmos y explica cuándo no hay solución. La tercera describe una sinusoide $A\,\mathrm{sen}(\omega t+\varphi)+k$: amplitud, periodo, frecuencia y valor inicial.
#
# **Convención.** En este libro $\ln$ es el logaritmo natural (base $e$) y $\log_{10}$ el decimal. En Python, `math.log(x)` y `np.log(x)` son el logaritmo **natural**; el decimal es `math.log10(x)`. Además, `math.sin`, `math.cos` y `np.sin` reciben **radianes**.

# %%
def trig_en(angulo, unidad):
    """sen, cos y tan de un ángulo. La unidad es obligatoria: 'grados' o 'radianes'."""
    if unidad not in ("grados", "radianes"):
        raise ValueError("La unidad debe ser 'grados' o 'radianes'.")
    t = math.radians(angulo) if unidad == "grados" else angulo
    s, c = math.sin(t), math.cos(t)
    s = 0.0 if abs(s) < 1e-12 else s          # lo menor que 1e-12 es ruido numérico, no un valor
    c = 0.0 if abs(c) < 1e-12 else c
    tg = s / c if c != 0 else None            # tan no existe donde cos = 0
    return s, c, tg


def calculadora_angulos(angulo, unidad, n_cifras):
    try:
        s, c, tg = trig_en(angulo, unidad)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    rad = math.radians(angulo) if unidad == "grados" else angulo
    print(f"ángulo = {cifras(angulo, n_cifras)} {unidad}  =  {cifras(rad, n_cifras)} rad  =  {cifras(math.degrees(rad), n_cifras)}°")
    print("sen =", cifras(s, n_cifras), "     cos =", cifras(c, n_cifras),
          "     tan =", cifras(tg, n_cifras) if tg is not None else "no existe (cos = 0)")
    print("sen² + cos² =", cifras(s**2 + c**2, n_cifras))


widgets.interact(calculadora_angulos,
    angulo=widgets.FloatText(value=30.0, description="ángulo"),
    unidad=widgets.Dropdown(options=["grados", "radianes"], value="grados", description="unidad"),
    n_cifras=widgets.IntSlider(value=5, min=1, max=8, description="cifras sig."));

# %%
def resolver_exponencial(base, c):
    """Resuelve base**x = c. Devuelve x o lanza ValueError con la razón."""
    if base <= 0 or base == 1:
        raise ValueError("La base debe cumplir base > 0 y base ≠ 1.")
    if c <= 0:
        raise ValueError("base**x siempre es positiva: con c ≤ 0 no hay solución real.")
    return math.log(c) / math.log(base)


def calculadora_exponencial(base, c, n_cifras):
    try:
        xs = resolver_exponencial(base, c)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    b_, c_ = cifras(base, n_cifras), cifras(c, n_cifras)
    print(f"{b_}^x = {c_}   →   x = ln({c_}) / ln({b_}) = {cifras(xs, n_cifras)}")
    print("comprobación: base^x =", cifras(base**xs, n_cifras))
    xx = np.linspace(xs - 3, xs + 3, 300)
    plt.figure(figsize=(4.5, 3)); plt.plot(xx, base**xx, "b-")
    plt.axhline(c, color="gray", ls="--"); plt.plot([xs], [c], "rs")
    plt.axhline(0, color="k", lw=0.8); plt.grid(True); plt.xlabel("x"); plt.ylabel("base^x"); plt.show()


widgets.interact(calculadora_exponencial,
    base=widgets.FloatText(value=5.0, description="base"), c=widgets.FloatText(value=7.0, description="c"),
    n_cifras=widgets.IntSlider(value=5, min=1, max=8, description="cifras sig."));

# %%
def sinusoide(A, w, fi, k=0.0):
    """Parámetros de v(t) = A sen(w t + fi) + k, con w en rad/s y fi en rad."""
    if A == 0:
        raise ValueError("A = 0: la señal es constante y no tiene amplitud.")
    if w == 0:
        raise ValueError("ω = 0: la señal es constante y no tiene periodo.")
    return {"amplitud": abs(A), "periodo": 2*math.pi/abs(w), "frecuencia": abs(w)/(2*math.pi),
            "v(0)": A*math.sin(fi) + k, "máximo": k + abs(A), "mínimo": k - abs(A)}


def calculadora_sinusoide(A, w, fi, k, n_cifras):
    try:
        r = sinusoide(A, w, fi, k)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"amplitud = {cifras(r['amplitud'], n_cifras)}      periodo T = {cifras(r['periodo'], n_cifras)} s      "
          f"frecuencia f = {cifras(r['frecuencia'], n_cifras)} Hz")
    print(f"v(0) = {cifras(r['v(0)'], n_cifras)}      máximo = {cifras(r['máximo'], n_cifras)}      mínimo = {cifras(r['mínimo'], n_cifras)}")
    t = np.linspace(0, 2*r["periodo"], 600)
    plt.figure(figsize=(5, 3)); plt.plot(t, A*np.sin(w*t + fi) + k, "b-")
    plt.axhline(k, color="gray", ls="--"); plt.axhline(0, color="k", lw=0.8); plt.grid(True)
    plt.xlabel("t (s)"); plt.ylabel("v(t)"); plt.show()


widgets.interact(calculadora_sinusoide,
    A=widgets.FloatText(value=8.0, description="A"), w=widgets.FloatText(value=100*math.pi, description="ω (rad/s)"),
    fi=widgets.FloatText(value=math.pi/4, description="φ (rad)"), k=widgets.FloatText(value=0.0, description="k"),
    n_cifras=widgets.IntSlider(value=5, min=1, max=8, description="cifras sig."));

# %%
def cerca(u, v, tol=1e-9):
    return abs(u - v) < tol

PRUEBAS_6 = [
    ("sen 30° = 0.5 y cos 60° = 0.5",
     lambda: cerca(trig_en(30, "grados")[0], 0.5) and cerca(trig_en(60, "grados")[1], 0.5)),
    ("30 radianes no es 30 grados: sen(30 rad) = -0.988",
     lambda: cerca(trig_en(30, "radianes")[0], -0.9880316, 1e-6)),
    ("tan 90° no existe", lambda: trig_en(90, "grados")[2] is None),
    ("la unidad es obligatoria", lambda: _rechaza(lambda: trig_en(30, "gradianes"))),
    ("5^x = 7 da x = 1.2091", lambda: round(resolver_exponencial(5, 7), 4) == 1.2091),
    ("2^x = 40 da x = log2(40) = 5.3219", lambda: round(resolver_exponencial(2, 40), 4) == 5.3219),
    ("c ≤ 0 se rechaza con mensaje", lambda: _rechaza(lambda: resolver_exponencial(2, -3))),
    ("base 1 se rechaza con mensaje", lambda: _rechaza(lambda: resolver_exponencial(1, 5))),
    ("8 sen(100πt + π/4): T = 0.02 s, f = 50 Hz, v(0) = 5.657",
     lambda: (lambda r: cerca(r["periodo"], 0.02) and cerca(r["frecuencia"], 50.0)
              and round(r["v(0)"], 3) == 5.657)(sinusoide(8, 100*math.pi, math.pi/4))),
    ("dominio de ln(4 - x²) es (-2, 2)",
     lambda: continuous_domain(parsear("log(4 - x**2)"), x, sp.S.Reals) == sp.Interval.open(-2, 2)),
    ("log2(40) - log2(5) = 3", lambda: cerca(math.log2(40) - math.log2(5), 3.0)),
]

for nombre, prueba in PRUEBAS_6:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 6).** Calcula a mano $\log_2 40-\log_2 5$ y $\mathrm{sen}\,150^\circ$ (usa el ángulo de referencia). Escribe tus resultados y compáralos con la calculadora. Anota: ¿usaste la propiedad del cociente? ¿de qué cuadrante salió el signo?

# %%
mi_log = None          # por ejemplo: 3
mi_seno = None         # por ejemplo: 0.5

esperado_log = math.log2(40) - math.log2(5)
esperado_seno = trig_en(150, "grados")[0]
print("La calculadora da: log2(40) - log2(5) =", round(esperado_log, 6), "| sen 150° =", round(esperado_seno, 6))
for nombre, mio, ref in (("log2(40) - log2(5)", mi_log, esperado_log), ("sen 150°", mi_seno, esperado_seno)):
    if mio is None:
        print(f"{nombre}: falta tu cálculo a mano.")
    else:
        print(f"{nombre}:", "coincide" if cerca(mio, ref, 1e-6) else "NO coincide, revisa tu paso")
