# ID: INT-U6-NB08
# Notebook: int/u6_aplicaciones.ipynb · sección 6.8 fuerza hidrostática y caudal
# Repositorio: int/u6_aplicaciones/08_fluidos.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6.8 Fuerza hidrostática y caudal
#
# La presión del agua a la profundidad $h$ es $p=\rho g h$. Sobre una placa vertical cuyo ancho a la profundidad $h$ es $w(h)$, entre $h_1$ y $h_2$:
#
# $$F=\int_{h_1}^{h_2}\rho g\,h\,w(h)\,dh,\qquad h_p=\frac1F\int_{h_1}^{h_2}h\cdot\rho g\,h\,w(h)\,dh\quad(\text{profundidad del centro de presión}).$$
#
# El centro de presión queda más abajo que el centroide de la placa, porque abajo empuja más. En un tubo circular de radio $R$ con perfil de velocidad $v(r)$, cada anillo de radio $r$ y ancho $dr$ tiene área $2\pi r\,dr$:
#
# $$Q=\int_0^Rv(r)\,2\pi r\,dr\quad(\text{m}^3/\text{s}),\qquad \bar v=\frac{Q}{\pi R^2}.$$
#
# El ejemplo exacto es una compuerta rectangular de 2 m de ancho entre 1 m y 2.5 m de profundidad, con $\rho=1000$ kg/m³ y $g=9.81$ m/s².

# %%
rho_g = 1000 * sp.Rational("9.81")
F_ej = sp.integrate(rho_g * x * 2, (x, 1, sp.Rational(5, 2)))
hp_ej = sp.integrate(x * rho_g * x * 2, (x, 1, sp.Rational(5, 2))) / F_ej
print("F =", F_ej, "N ≈", cifras(F_ej, 6), "N   h_p =", hp_ej, "m ≈", cifras(hp_ej, 7), "m")

# %% [markdown]
# **Calculadora de la fuerza sobre una placa vertical.** Escribe el ancho $w(h)$ en m como función de la profundidad $h$ (constante para un rectángulo, `h` para un triángulo con el vértice arriba), las profundidades del borde superior y del inferior y la densidad del líquido. La calculadora da la fuerza, el centroide y el centro de presión, y dibuja la placa con la profundidad hacia abajo.

# %%
def fuerza_placa(w, h1, h2, rho=1000, g=9.81):
    """Fuerza hidrostática (N) sobre una placa vertical de ancho w(h) entre las profundidades h1 < h2 (m), y la
    profundidad h_p (m) del centro de presión. También devuelve el área A y la profundidad del centroide h_c."""
    w = sp.sympify(w)
    h1, h2 = exacto(h1, "h1"), exacto(h2, "h2")
    rho, g = exacto(rho, "ρ"), exacto(g, "g")
    if h1 < 0:
        raise ValueError("h1 es una profundidad bajo la superficie: no puede ser negativa.")
    if not h1 < h2:
        raise ValueError("el borde inferior (h2) debe estar más hondo que el superior (h1).")
    if rho <= 0 or g <= 0:
        raise ValueError("ρ y g deben ser positivas.")
    revisa_intervalo(w, h1, h2, "w", var="h")
    anchos = numerica(w)(malla_interior(h1, h2))
    if np.any(anchos < -1e-12):
        raise ValueError("el ancho w(h) sale negativo en parte del intervalo: revisa la fórmula.")
    integ = Integrador()
    def _int(f):
        return integ(sp.expand(f), h1, h2)
    A = _int(w)
    if float(A) <= 0:
        raise ValueError("la placa tiene área cero.")
    F = simplifica(rho * g * _int(x * w))
    hp = simplifica(rho * g * _int(x**2 * w) / F)
    hc = simplifica(_int(x * w) / A)
    return {"F": F, "h_p": hp, "A": A, "h_c": hc, "metodo": integ.metodo}


def calculadora_placa(w_txt, h1, h2, rho):
    try:
        w = parsear(w_txt)
        r = fuerza_placa(w, h1, h2, rho)
    except ERRORES as err:
        print("Revisa la entrada:", explica(err)); return
    print(f"({CIFRAS} cifras significativas, redondeadas; g = 9.81 m/s²; {COMO[r['metodo']]})")
    print(f"F   = {muestra(r['F'])} N = {cifras(float(r['F']) / 1000, CIFRAS)} kN")
    print(f"A   = {muestra(r['A'])} m²    centroide a h_c = {muestra(r['h_c'])} m")
    print(f"centro de presión a h_p = {muestra(r['h_p'])} m de profundidad")
    try:
        hs = np.linspace(float(h1), float(h2), 300)
        ws = para_graficar(w, hs)
        wmax = float(np.nanmax(ws))
        fig, ax = plt.subplots(figsize=(4.8, 3.4))
        ax.fill_betweenx(hs, -ws / 2, ws / 2, facecolor="0.85", edgecolor="black")
        ax.axhline(0, color="black", lw=1.2); rotula(ax, -1.5 * wmax, 0, "superficie", va="bottom", xytext=(0, 3))
        for prof, estilo, texto, lado in ((float(r["h_c"]), ":", "centroide", -1),
                                          (float(r["h_p"]), "--", "centro de presión", 1)):
            ax.axhline(prof, color="black", ls=estilo, lw=1)
            rotula(ax, lado * wmax / 2, prof, f"{texto}, {cifras(prof, CIFRAS)} m", xytext=(lado * 6, 3 * lado),
                   ha="left" if lado > 0 else "right", va="top" if lado > 0 else "bottom")
        ax.set_ylim(float(h2) * 1.08, -0.1 * float(h2)); ax.set_xlim(-1.6 * wmax, 1.6 * wmax)
        ax.set_xlabel("ancho (m)"); ax.set_ylabel("profundidad h (m)"); plt.show()
    except ERRORES:
        print("(no pude dibujar la gráfica con estos datos)")


widgets.interact(calculadora_placa,
    w_txt=widgets.Text(value="2", description="w(h) =", continuous_update=False),
    h1=widgets.FloatText(value=1, description="h1 (m)"),
    h2=widgets.FloatText(value=2.5, description="h2 (m)"),
    rho=widgets.FloatText(value=1000, description="ρ (kg/m³)"));

# %% [markdown]
# **Calculadora del caudal en un tubo.** Escribe el perfil de velocidad $v(r)$ en m/s con la variable `r` y el radio $R$ del tubo en m. El perfil del ejemplo es el laminar, $v=v_{\max}\left(1-r^2/R^2\right)$ con $v_{\max}=2$ m/s; prueba también uno turbulento, `(1 - r/0.05)^(1/7)`. La calculadora da el caudal y la velocidad media y la compara con la máxima.

# %%
def caudal(v, R):
    """Caudal Q = ∫_0^R v(r) 2πr dr (m³/s) y velocidad media Q/(πR²) (m/s) en un tubo circular de radio R (m)."""
    v = sp.sympify(v)
    R = exacto(R, "R")
    if R <= 0:
        raise ValueError("el radio R debe ser positivo.")
    revisa_intervalo(v, 0, R, "v", var="r")
    rs = np.linspace(0, float(R), 403)
    vs = numerica(v)(rs)
    if not np.all(np.isfinite(vs)):
        raise ValueError("la velocidad debe ser finita en todo el tubo, también en el centro (r = 0) y en la pared.")
    if np.any(vs < -1e-12):
        raise ValueError(f"v(r) sale negativa cerca de r ≈ {cifras(rs[int(np.argmin(vs))], 3)} m: revisa que el radio "
                         "que usaste en la fórmula coincida con R.")
    Q, como = integra(sp.expand(2 * sp.pi * x * v), 0, R)
    Q = simplifica(Q)
    if float(Q) <= 0:
        raise ValueError("el caudal es cero: con v(r) = 0 en todo el tubo no hay flujo.")
    return {"Q": Q, "v_media": simplifica(Q / (sp.pi * R**2)), "metodo": como}


def calculadora_caudal(v_txt, R):
    try:
        v = parsear(v_txt)
        r = caudal(v, R)
        rs = np.linspace(0, R, 300)
        vs = numerica(v)(rs)
    except ERRORES as err:
        print("Revisa la entrada:", explica(err)); return
    print(f"({CIFRAS} cifras significativas, redondeadas; {COMO[r['metodo']]})")
    print(f"Q = {muestra(r['Q'])} m³/s = {cifras(1000 * float(r['Q']), CIFRAS)} L/s")
    print(f"velocidad media = {muestra(r['v_media'])} m/s   máxima = {cifras(float(np.max(vs)), CIFRAS)} m/s   "
          f"media/máxima = {cifras(float(r['v_media']) / float(np.max(vs)), 4)}  (cociente con 4 cifras)")
    try:
        fig, ax = plt.subplots(figsize=(4.8, 3.2))
        rr = np.concatenate((-rs[::-1], rs)); vv = np.concatenate((vs[::-1], vs))
        ax.plot(vv, rr, color="black"); ax.fill_betweenx(rr, 0, vv, color="0.88")
        vm = float(r["v_media"])
        ax.axvline(vm, color="black", ls="--", lw=1); rotula(ax, vm, -0.8 * R, f"media {cifras(vm, CIFRAS)} m/s")
        ax.axhline(R, color="black", lw=2); ax.axhline(-R, color="black", lw=2)
        ax.set_xlabel("v (m/s)"); ax.set_ylabel("r (m)"); ax.set_title("perfil de velocidad (pared en ±R)", fontsize=9)
        plt.show()
    except ERRORES:
        print("(no pude dibujar la gráfica con estos datos)")


widgets.interact(calculadora_caudal,
    v_txt=widgets.Text(value="2*(1 - r^2/0.05^2)", description="v(r) =", continuous_update=False),
    R=widgets.FloatText(value=0.05, description="R (m)"));

# %%
# Casos de prueba de la sección 6.8 (resultado conocido)
PRUEBAS_8 = [
    ("compuerta de 2 m de ancho entre 1 y 2.5 m: F = 51502.5 N y h_p = 1.857143 m",
     lambda: (lambda r: r["F"] == sp.Rational("51502.5") and cifras(r["h_p"], 7) == "1.857143")(fuerza_placa(2, 1, 2.5))),
    ("placa circular de radio 1 m con el centro a 1 m: F = ρ g π = 9810π N y h_p = 1 + (π/4)/π = 1.25 m",
     lambda: (lambda r: sp.simplify(r["F"] - 9810 * sp.pi) == 0 and r["h_p"] == sp.Rational(5, 4))(
         fuerza_placa(2 * sp.sqrt(1 - (x - 1)**2), 0, 2))),
    ("placa cuadrada de 2 m con el borde en la superficie: 39240 N y h_p = 4/3 m",
     lambda: (lambda r: r["F"] == 39240 and r["h_p"] == sp.Rational(4, 3))(fuerza_placa(2, 0, 2))),
    ("placa triangular w = h de 0 a 3 m: 88290 N", lambda: fuerza_placa(x, 0, 3)["F"] == 88290),
    ("borde inferior arriba del superior se rechaza", lambda: _rechaza(lambda: fuerza_placa(1, 2, 1))),
    ("perfil laminar 2(1 - r²/R²), R = 0.05 m: velocidad media 1 m/s",
     lambda: caudal(2 * (1 - x**2 / sp.Rational(1, 20)**2), 0.05)["v_media"] == 1),
    ("perfil (1 - r/R)^(1/7), R = 0.1 m: velocidad media 49/60 de la máxima",
     lambda: caudal((1 - x / sp.Rational(1, 10))**sp.Rational(1, 7), 0.1)["v_media"] == sp.Rational(49, 60)),
    ("perfil uniforme v = 3 m/s: velocidad media 3 m/s", lambda: caudal(3, 0.2)["v_media"] == 3),
    ("radio negativo se rechaza", lambda: _rechaza(lambda: caudal(3, -0.1))),
    ("perfil con velocidad negativa (radio de la fórmula distinto de R) se rechaza",
     lambda: _rechaza(lambda: caudal(2 * (1 - x**2 / sp.Rational(1, 20)**2), 0.1))),
]
for nombre, prueba in PRUEBAS_8:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 6.8).** Una placa cuadrada de 1 m de lado está vertical en agua, con el borde superior en la superficie. Calcula a mano la fuerza del agua sobre una cara ($\rho=1000$ kg/m³, $g=9.81$ m/s²) y escribe el número en N.

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
        ref = fuerza_placa(1, 0, 1)["F"]
        print("La calculadora da: F =", muestra(ref), "N")
        print("coinciden a 3 cifras" if coincide(mio, ref) else
              "NO coinciden: F = ∫_0^1 ρ g h · 1 dh = ρ g/2; la presión crece con la profundidad, no es ρ g · 1 en toda la placa")
