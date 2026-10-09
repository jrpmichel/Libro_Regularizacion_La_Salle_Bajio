# ID: INT-U5-NB04
# Notebook: int/u5_tecnicas.ipynb · sección 5.4 potencias trigonométricas y sustitución trigonométrica
# Repositorio: int/u5_tecnicas/04_trigonometricas.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 4. Potencias trigonométricas y sustitución trigonométrica
#
# Identidades de reducción: $\;\sin^2x=\frac{1-\cos2x}{2},\quad\cos^2x=\frac{1+\cos2x}{2},\quad\sin x\cos x=\frac12\sin2x.$
#
# Para $\int\sin^mx\cos^nx\,dx$: si $m$ es impar se usa $u=\cos x$; si $n$ es impar, $u=\sin x$; si los dos son pares, se bajan las potencias.
#
# Sustitución trigonométrica ($a>0$): $\sqrt{a^2-x^2}$ con $x=a\sin\theta$; $\sqrt{a^2+x^2}$ con $x=a\tan\theta$; $\sqrt{x^2-a^2}$ con $x=a/\cos\theta$.

# %%
print("∫_0^{2π} sen² x dx =", sp.integrate(sp.sin(x)**2, (x, 0, 2 * sp.pi)), " (la mitad del rectángulo 2π × 1)")
a_ = sp.symbols("a", positive=True)
print("∫_0^a √(a² - x²) dx =", sp.integrate(sp.sqrt(a_**2 - x**2), (x, 0, a_)), " (cuarto de círculo)")

# %% [markdown]
# **Calculadora de potencias.** Escribe los exponentes $m$ y $n$ de $\sin^mx\cos^nx$. Te dice qué caso aplica, muestra el integrando reescrito y el resultado.

# %%
from sympy.simplify.fu import TR8 as _TR8      # productos y potencias de seno y coseno -> sumas


def potencias_trig(m, n):
    """(caso, integrando reescrito, F) para ∫ sen^m x cos^n x dx, con m y n enteros no negativos.
    En los casos (a) y (b) el integrando reescrito está en u; en el (c), en x, ya sin potencias (TR8)."""
    if int(m) != m or int(n) != n or m < 0 or n < 0:
        raise ValueError("m y n deben ser enteros no negativos.")
    m, n = int(m), int(n)
    if m % 2 == 1:
        caso = "(a) m impar: separa un sen x, cambia sen² = 1 - cos² y usa u = cos x (du = -sen x dx)"
        reescrita = sp.expand(-(1 - U**2)**((m - 1) // 2) * U**n)
        F = sp.integrate(reescrita, U).subs(U, sp.cos(x))
    elif n % 2 == 1:
        caso = "(b) n impar: separa un cos x, cambia cos² = 1 - sen² y usa u = sen x (du = cos x dx)"
        reescrita = sp.expand(U**m * (1 - U**2)**((n - 1) // 2))
        F = sp.integrate(reescrita, U).subs(U, sp.sin(x))
    else:
        caso = "(c) los dos pares: baja las potencias con las identidades de reducción"
        reescrita = sp.sin(x)**m * sp.cos(x)**n
        for _ in range(8):                     # TR8 repetido hasta que solo queden cosenos de múltiplos de x
            nueva = sp.expand(_TR8(sp.expand(reescrita)))
            if nueva == reescrita:
                break
            reescrita = nueva
        F = sp.integrate(reescrita, x)
    return caso, reescrita, sp.expand(F)


def calculadora_potencias(m, n):
    try:
        caso, reescrita, F = potencias_trig(m, n)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    f = sp.sin(x)**m * sp.cos(x)**n
    print("∫", f, "dx")
    print("caso:", caso)
    print("integrando reescrito:", reescrita, "(en x)" if m % 2 == 0 and n % 2 == 0 else "(en u, con du ya incluido)")
    print("resultado:", mostrar_antiderivada(F))
    print("derivada del resultado = integrando:", "sí" if compara_con_sympy(F, f)[0] else "NO")


widgets.interact(calculadora_potencias,
    m=widgets.IntSlider(value=3, min=0, max=6, description="m (sen)", continuous_update=False),
    n=widgets.IntSlider(value=0, min=0, max=6, description="n (cos)", continuous_update=False));

# %% [markdown]
# **Calculadora de sustitución trigonométrica.** Escribe una función con una raíz $\sqrt{a^2-x^2}$, $\sqrt{a^2+x^2}$ o $\sqrt{x^2-a^2}$ (también sirven potencias como $(x^2+9)^{3/2}$). La calculadora reconoce el caso, propone la sustitución, dibuja el triángulo, integra en $\theta$ y regresa a $x$ con las razones del triángulo. En el caso $\sqrt{x^2-a^2}$ trabaja con $x>a$.

# %%
TH = sp.symbols("theta", real=True)
_S, _K = sp.sin(TH), sp.cos(TH)


def _raiz_cuadratica(f):
    """La potencia (±x² + c)^(k/2) de f, con su signo de x² y su constante c; ValueError si hay dos raíces distintas."""
    bases = []
    for p in f.atoms(sp.Pow):
        if p.exp.is_Rational and p.exp.q == 2 and p.base.is_polynomial(x):
            b = sp.expand(p.base)
            if b not in bases:
                bases.append(b)
    if len(bases) > 1:
        raise ValueError("f tiene dos raíces distintas; la sustitución trigonométrica trabaja con una sola.")
    for b in bases:
        coef = sp.Poly(b, x).all_coeffs()
        if len(coef) == 3 and coef[1] == 0 and coef[0] in (1, -1) and coef[2] != 0:
            return b, coef[0], coef[2]
    return None


def _razones(caso, a):
    """Triángulo de referencia: sen θ, cos θ y tan θ escritos con x, y θ como función de x."""
    if caso == "a² - x²":       # x = a sen θ: opuesto x, hipotenusa a, adyacente √(a² - x²)
        return {_S: x / a, _K: sp.sqrt(a**2 - x**2) / a, sp.tan(TH): x / sp.sqrt(a**2 - x**2)}, sp.asin(x / a)
    if caso == "a² + x²":       # x = a tan θ: opuesto x, adyacente a, hipotenusa √(a² + x²)
        return {sp.tan(TH): x / a, _K: a / sp.sqrt(a**2 + x**2), _S: x / sp.sqrt(a**2 + x**2)}, sp.atan(x / a)
    # x = a/cos θ (x > a): hipotenusa x, adyacente a, opuesto √(x² - a²)
    return {_K: a / x, _S: sp.sqrt(x**2 - a**2) / x, sp.tan(TH): sp.sqrt(x**2 - a**2) / a}, sp.acos(a / x)


def _forma_estandar(H):
    """Reescribe H(θ) solo con sen θ, cos θ y θ. sympy da ∫ sec θ dθ como ½ln(sen θ + 1) - ½ln(sen θ - 1);
    como cos²θ = (1 - sen θ)(1 + sen θ), eso es ln|sec θ + tan θ|, y con el triángulo queda real.
    Igual con ∫ csc θ dθ = ln|csc θ - cot θ|."""
    H = H.replace(lambda e: isinstance(e, sp.tan) and e.args[0] == TH / 2, lambda e: _S / (1 + _K))
    H = sp.expand_trig(sp.expand(H))
    for menos, mas, otra in ((_S, 1 + _S, _K), (_K, 1 + _K, _S)):
        H = H.replace(lambda e: isinstance(e, sp.log) and sp.expand(e.args[0] ** 2 - (menos - 1) ** 2) == 0,
                      lambda e: 2 * sp.log(otra) - sp.log(mas))     # |sen θ - 1| = cos²θ/(1 + sen θ)
    return sp.expand(H)


def _sin_constantes(F):
    """Quita los sumandos que no dependen de x (como -ln 3): F sigue siendo antiderivada."""
    F = sp.expand(sp.expand_log(F, force=True))
    return sp.Add(*[t for t in sp.Add.make_args(F) if t.has(x)])


def _ordena(F):
    """Simplifica por separado la parte algebraica y cada logaritmo, para que no se mezclen en una fracción."""
    F = sp.expand(F)
    logs = [t for t in sp.Add.make_args(F) if t.has(sp.log)]
    resto = [t for t in sp.Add.make_args(F) if not t.has(sp.log)]
    return sp.simplify(sp.Add(*resto)) + sp.Add(*[sp.simplify(t) for t in logs])


def sustitucion_trig(f):
    """(caso, a, x(θ), raíz en θ, integrando en θ, H(θ), F(x)); ValueError si no hay una raíz de esas formas
    o si el resultado no se puede regresar a x con el triángulo."""
    r = _raiz_cuadratica(f)
    if r is None:
        raise ValueError("no encontré √(a² - x²), √(a² + x²) ni √(x² - a²) (con x² de coeficiente 1 y sin término en x).")
    base, signo, c = r
    if signo == -1 and c > 0:
        a = sp.sqrt(c); caso = "a² - x²"; xt, raiz = a * _S, a * _K
    elif signo == 1 and c > 0:
        a = sp.sqrt(c); caso = "a² + x²"; xt, raiz = a * sp.tan(TH), a / _K
    elif signo == 1 and c < 0:
        a = sp.sqrt(-c); caso = "x² - a²"; xt, raiz = a / _K, a * sp.tan(TH)
    else:
        raise ValueError("la expresión bajo la raíz es negativa para todo x.")
    g = f.replace(lambda e: isinstance(e, sp.Pow) and sp.expand(e.base) == base and e.exp.is_Rational and e.exp.q == 2,
                  lambda e: raiz**(2 * e.exp))
    g = sp.simplify(g.subs(x, xt) * sp.diff(xt, TH))
    if g.has(x):
        raise ValueError("después de sustituir queda una x: revisa que la raíz sea la única parte con x² + c.")
    H = con_limite(lambda: sp.integrate(g, TH), 15)
    if H.has(sp.Integral):
        raise ValueError(f"sympy no encontró ∫ {g} dθ con funciones elementales.")
    H = _forma_estandar(H)
    razones, inversa = _razones(caso, a)
    F = H.subs(razones).subs(TH, inversa)
    F = con_valor_absoluto(_ordena(_sin_constantes(sp.simplify(F))))
    puntos = _puntos_dominio(caso, a)
    if not comprueba_en(F, f, puntos):
        raise ValueError("el resultado no se pudo regresar a x con el triángulo (o no es real en el dominio).")
    return caso, a, xt, raiz, g, H, F


def _puntos_dominio(caso, a):
    """Puntos de prueba dentro del dominio (sin x = 0, donde algunos integrandos como √(a² + x²)/x no existen)."""
    a = float(a)
    return {"a² - x²": [-0.7 * a, -0.3 * a, 0.2 * a, 0.5 * a, 0.8 * a],
            "a² + x²": [-2 * a, -0.5 * a, 0.3 * a, a, 2.5 * a],
            "x² - a²": [1.2 * a, 1.5 * a, 2 * a, 3 * a, 4 * a]}[caso]


def comprueba_en(F, f, puntos):
    """True si F es real y F' = f en todos los puntos (puntos dentro del dominio de f)."""
    dF = derivada_legible(F)
    for p in puntos:
        try:
            if abs(valor_real(F, p)) > 1e12 or abs(valor_real(dF, p) - valor_real(f, p)) > 1e-9:
                return False
        except (ValueError, TypeError, ZeroDivisionError, OverflowError):
            return False
    return True


def dibuja_triangulo(caso, a):
    lados = {"a² - x²": ("x", "√(a² - x²)", "a"), "a² + x²": ("x", "a", "√(a² + x²)"),
             "x² - a²": ("√(x² - a²)", "a", "x")}[caso]           # (opuesto, adyacente, hipotenusa)
    fig, ax = plt.subplots(figsize=(3, 2.2))
    ax.plot([0, 2, 2, 0], [0, 0, 1.3, 0], color="black", lw=1.5)
    ax.text(2.06, 0.6, lados[0]); ax.text(0.85, -0.18, lados[1]); ax.text(0.7, 0.78, lados[2])
    ax.text(0.35, 0.06, "θ"); ax.set_axis_off(); ax.set_title(f"caso {caso}, a = {a}", fontsize=9)
    plt.show()


def calculadora_sustitucion_trig(f_txt, dibujar):
    try:
        f = parsear(f_txt)
        caso, a, xt, raiz, g, H, F = sustitucion_trig(f)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"caso {caso} con a = {a}:  x = {xt},  la raíz vale {raiz}" + ("   (para x > a)" if caso == "x² - a²" else ""))
    print("integrando en θ:", g)
    print("integral en θ:", H, "+ C")
    print("con el triángulo, de regreso a x:", mostrar_antiderivada(F))
    print("derivada del resultado = f(x) en el dominio: sí (comprobado en 5 puntos)")
    if dibujar:
        dibuja_triangulo(caso, a)


widgets.interact(calculadora_sustitucion_trig,
    f_txt=widgets.Text(value="x^2/sqrt(4-x^2)", description="f(x) =", continuous_update=False),
    dibujar=widgets.Checkbox(value=True, description="dibujar el triángulo"));

# %% [markdown]
# **Calculadora de valor eficaz.** $V_{\mathrm{rms}}=\sqrt{\frac1T\int_0^Tv^2\,dt}$. Escribe $v(t)$ y el periodo $T$ en segundos (puedes escribir `1/60` o `2*pi`). Con *media onda*, la señal vale $v(t)$ en la primera mitad del periodo y $0$ en la segunda.

# %%
def valor_eficaz(v, T, media_onda=False):
    """√((1/T) ∫_0^T v² dt); con media_onda, v solo cuenta en [0, T/2]. T puede ser texto ("1/60", "2*pi")."""
    Tn = parsear(T) if isinstance(T, str) else sp.nsimplify(T, [sp.pi])
    if Tn.free_symbols or not Tn.is_real or not Tn > 0:
        raise ValueError("el periodo debe ser un número positivo.")
    fin = Tn / 2 if media_onda else Tn
    for k in (0.13, 0.37, 0.71, 0.93):                  # v debe ser real en el intervalo
        try:
            valor_real(v, k * fin)
        except ValueError:
            raise ValueError(f"v no tiene valor real en t = {cifras(float(k * fin), 3)} s: revisa la función.") from None
    I = con_limite(lambda: sp.integrate(v**2, (x, 0, fin)), 15)
    especiales = [a for a in I.atoms(sp.Function) if not isinstance(a, _PERMITIDAS)]
    if I.has(sp.Integral) or especiales:                # sin antiderivada elemental: cuadratura
        valor, err = quad(sp.lambdify(x, v**2, "math"), 0, float(fin), limit=200)
        if not math.isfinite(valor) or err > 1e-8 * max(1.0, abs(valor)):
            raise ValueError("la integral de v² no converge en el periodo: v no tiene valor eficaz finito.")
        return sp.Float(math.sqrt(valor / float(Tn)), 15)
    if not I.is_finite:
        raise ValueError("la integral de v² diverge en el periodo: v no tiene valor eficaz finito.")
    return sp.sqrt(I / Tn)


def calculadora_eficaz(v_txt, T_txt, media_onda, n_cifras):
    try:
        v = parsear(v_txt)
        r = valor_eficaz(v, T_txt, media_onda)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    exacto = "" if isinstance(r, sp.Float) else f" = {sp.nsimplify(sp.simplify(r))}"
    print(f"V_rms{exacto} ≈ {cifras(r, n_cifras)}   ({n_cifras} cifras significativas, redondeado; ángulos en radianes)")
    if isinstance(r, sp.Float):
        print("(v² no tiene antiderivada elemental: el valor se obtuvo con una suma numérica)")


widgets.interact(calculadora_eficaz,
    v_txt=widgets.Text(value="179.6*sin(120*pi*t)", description="v(t) =", continuous_update=False),
    T_txt=widgets.Text(value="1/60", description="T (s) =", continuous_update=False),
    media_onda=widgets.Checkbox(value=False, description="media onda"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig.", continuous_update=False));

# %%
# Casos de prueba de la sección 4 (resultado conocido)
PRUEBAS_4 = [
    ("sen³: caso (a), -cos x + cos³x/3", lambda: (lambda r: r[0].startswith("(a)") and
     sp.simplify(r[2] - (-sp.cos(x) + sp.cos(x)**3 / 3)) == 0)(potencias_trig(3, 0))),
    ("sen² cos³: caso (b), F' = f", lambda: (lambda r: r[0].startswith("(b)") and
     compara_con_sympy(r[2], sp.sin(x)**2 * sp.cos(x)**3)[0])(potencias_trig(2, 3))),
    ("sen² cos²: caso (c), 1/8 - cos 4x/8", lambda: (lambda r: r[0].startswith("(c)") and
     sp.simplify(r[1] - (sp.Rational(1, 8) - sp.cos(4 * x) / 8)) == 0)(potencias_trig(2, 2))),
    ("√(4 - x²): caso a² - x², x√(4 - x²)/2 + 2 arcsen(x/2)", lambda: (lambda r: r[0] == "a² - x²" and
     misma_funcion(r[6], x * sp.sqrt(4 - x**2) / 2 + 2 * sp.asin(x / 2)))(sustitucion_trig(sp.sqrt(4 - x**2)))),
    ("(x² + 9)^(-3/2): caso a² + x², x/(9√(x² + 9))", lambda: (lambda r: r[0] == "a² + x²" and
     misma_funcion(r[6], x / (9 * sp.sqrt(x**2 + 9))))(sustitucion_trig((x**2 + 9)**sp.Rational(-3, 2)))),
    ("1/√(x² + 4): resultado real, ln|x + √(x² + 4)|", lambda: (lambda r: not r[6].has(sp.I) and
     comprueba_en(r[6], 1 / sp.sqrt(x**2 + 4), [-3, -0.5, 0.4, 2, 7]))(sustitucion_trig(1 / sp.sqrt(x**2 + 4)))),
    ("valor eficaz de sen(120πt) con T = 1/60: 1/√2", lambda: valor_eficaz(sp.sin(120 * sp.pi * x), "1/60") == 1 / sp.sqrt(2)),
    ("valor eficaz de la rampa v = t con T = 1: 1/√3", lambda: valor_eficaz(x, "1") == 1 / sp.sqrt(3)),
    ("valor eficaz de v = 1 con media onda: 1/√2", lambda: valor_eficaz(sp.Integer(1), "2", True) == 1 / sp.sqrt(2)),
]
for nombre, prueba in PRUEBAS_4:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 4).** Calcula a mano $\int_0^{\pi}\sin^2x\,dx$ y escribe el número.

# %%
mi_valor = None     # escribe un número (también sirve "pi/2")

if mi_valor is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    ref = sp.integrate(sp.sin(x)**2, (x, 0, sp.pi))
    try:
        mio = a_numero(mi_valor)
    except ValueError as err:
        print("Revisa tu número:", err)
    else:
        print("La calculadora da:", ref, "≈", cifras(ref, 6))
        print("coinciden a 3 cifras" if coincide(mio, ref) else
              "NO coinciden: sen² x = (1 - cos 2x)/2, y cos 2x acumula cero en [0, π]")
