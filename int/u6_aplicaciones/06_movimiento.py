# ID: INT-U6-NB06
# Notebook: int/u6_aplicaciones.ipynb · sección 6.6 de la aceleración a la posición
# Repositorio: int/u6_aplicaciones/06_movimiento.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6.6 De la aceleración a la posición
#
# $$v(t)=v_0+\int_0^ta(\tau)\,d\tau,\qquad s(t)=s_0+\int_0^tv(\tau)\,d\tau.$$
#
# Con datos medidos cada $\Delta t$, las integrales se acumulan con la regla del trapecio: $v_{k+1}=v_k+\tfrac{\Delta t}{2}(a_k+a_{k+1})$ y lo mismo para $s$ a partir de $v$. Un acelerómetro con sesgo constante $b$ mide $a+b$; al integrar dos veces, el error en la posición crece como $e(t)=\tfrac12bt^2$, aunque $b$ sea pequeño. Por eso la navegación con acelerómetros necesita una referencia externa que la corrija.
#
# El ejemplo exacto parte del reposo en el origen con $a(t)=6t$ m/s².

# %%
t_s, tau_s = sp.symbols("t tau", nonnegative=True)
v_ej = sp.integrate(6 * tau_s, (tau_s, 0, t_s))
s_ej = sp.integrate(v_ej.subs(t_s, tau_s), (tau_s, 0, t_s))
print("v(t) =", v_ej, "m/s   s(t) =", s_ej, "m   s(2) =", s_ej.subs(t_s, 2), "m")

# %% [markdown]
# **Calculadora de movimiento.** Escribe $a(t)$ en m/s² (con `escalon(t - 2)` armas perfiles por tramos), la velocidad y la posición iniciales y el tiempo final. La calculadora integra dos veces con `sympy` y, aparte, simula un acelerómetro con sesgo $b$: integra con trapecios las muestras de $a+b$ y compara la posición medida con la verdadera. El ejemplo es el perfil de un vehículo que acelera 2 s, avanza 16 s a velocidad constante y frena 2 s.

# %%
_TAU = sp.Dummy("tau", real=True)


class SinFormaExacta(ValueError):
    """sympy no encontró v(t) o s(t) como fórmula; la integral puede existir y calcularse con muestras."""


def movimiento(a, v0=0, s0=0):
    """v(t) y s(t) exactas (expresiones de sympy en x, que hace de t) a partir de a(t), v(0) = v0 y s(0) = s0.
    ValueError si la integral desde t = 0 diverge; SinFormaExacta si sympy no encuentra la fórmula (en ese caso
    integra muestras de a(t) con integra_datos)."""
    a = sp.sympify(a)
    v0, s0 = exacto(v0, "v0"), exacto(s0, "s0")
    if a.has(sp.Heaviside):
        a = a.rewrite(sp.Piecewise)
    def _acumula(g, cual):
        try:
            r = con_limite(lambda: sp.integrate(g.subs(x, _TAU), (_TAU, 0, x)), 10)
        except (ValueError, NotImplementedError, TypeError, ArithmeticError):
            raise SinFormaExacta("sympy no encontró una fórmula exacta para v(t) y s(t); integra muestras de a(t) "
                                 "con integra_datos.") from None
        if r.has(sp.oo, -sp.oo, sp.zoo, sp.nan):
            raise ValueError(f"la integral de {cual} desde t = 0 diverge: {cual} crece sin límite cerca de t = 0, así "
                             "que el movimiento no está definido con esos datos.")
        if r.has(sp.Integral, sp.exp_polar, sp.hyper, sp.meijerg):
            raise SinFormaExacta("sympy no encontró una fórmula exacta para v(t) y s(t); integra muestras de a(t) "
                                 "con integra_datos.")
        return r if r.has(sp.Piecewise, sp.Min, sp.Max) else sp.expand(r)
    v = v0 + _acumula(a, "a(t)")
    s = s0 + _acumula(v, "v(t)")
    return v, s


def integra_datos(t, a, v0=0.0, s0=0.0):
    """Velocidad y posición acumuladas con la regla del trapecio a partir de muestras a(t_k) (arreglos de numpy)."""
    t, a = np.asarray(t, dtype=float), np.asarray(a, dtype=float)
    if t.ndim != 1 or t.shape != a.shape or len(t) < 2:
        raise ValueError("t y a deben ser listas del mismo tamaño, con al menos dos muestras.")
    if not (np.all(np.isfinite(t)) and np.all(np.isfinite(a))):
        raise ValueError("hay muestras que no son números finitos.")
    dt = np.diff(t)
    if np.any(dt <= 0):
        raise ValueError("los tiempos deben ir en aumento.")
    v = v0 + np.concatenate(([0.0], np.cumsum(dt * (a[1:] + a[:-1]) / 2)))
    s = s0 + np.concatenate(([0.0], np.cumsum(dt * (v[1:] + v[:-1]) / 2)))
    return v, s


def error_sesgo(b, t):
    """Error en la posición, en m, por un sesgo constante b (m/s²) del acelerómetro después de t segundos."""
    return b * t**2 / 2


def tiempo_para_error(b, e_max):
    """Tiempo (s) en el que el error ½ b t² llega a e_max (m)."""
    if b <= 0 or e_max <= 0:
        raise ValueError("el sesgo y el error máximo deben ser positivos.")
    return math.sqrt(2 * e_max / b)


def _por_tramos(e):
    return e.has(sp.Piecewise, sp.Min, sp.Max, sp.Heaviside)


def calculadora_movimiento(a_txt, v0, s0, t_fin, b):
    try:
        a = parsear(a_txt)
        if not t_fin > 0:
            raise ValueError("el tiempo final debe ser positivo.")
        tf = exacto(t_fin, "el tiempo final")
        revisa_intervalo(a, 0, tf, "a(t)", var="t")
        try:
            v, s = movimiento(a, v0, s0)
        except SinFormaExacta:
            v = s = None
        ts = np.linspace(0, t_fin, 2001)
        a_k = numerica(a)(ts)
    except ERRORES as err:
        print("Revisa la entrada:", explica(err)); return
    t_sim = sp.Symbol("t")
    print(f"({CIFRAS} cifras significativas, redondeadas; t en s, a en m/s², v en m/s, s en m)")
    if v is None:
        print("sympy no encontró v(t) y s(t) como fórmula; los valores salen de integrar 2001 muestras con trapecios")
    elif _por_tramos(v) or _por_tramos(s):
        print("v(t) y s(t) quedan por tramos; la gráfica muestra s(t)")
    else:
        print(f"v(t) = {v.subs(x, t_sim)}     s(t) = {s.subs(x, t_sim)}")
    if v is not None:
        print(f"en t = {punto(tf)} s:  v = {muestra(simplifica(v.subs(x, tf)))} m/s   "
              f"s = {muestra(simplifica(s.subs(x, tf)))} m")
    if not np.all(np.isfinite(a_k)):
        print("a(t) no es finita en t = 0: se omite la simulación del acelerómetro, que usa muestras desde t = 0.")
        return
    v_ver, s_ver = integra_datos(ts, a_k, v0, s0)
    _, s_med = integra_datos(ts, a_k + b, v0, s0)
    if v is None:
        print(f"en t = {punto(tf)} s:  v ≈ {cifras(v_ver[-1], CIFRAS)} m/s   s ≈ {cifras(s_ver[-1], CIFRAS)} m")
    print(f"acelerómetro con sesgo b = {b:g} m/s²: posición medida {cifras(s_med[-1], CIFRAS)} m, "
          f"error {cifras(s_med[-1] - s_ver[-1], CIFRAS)} m  (½ b t² = {cifras(error_sesgo(b, t_fin), CIFRAS)} m)")
    if b > 0:
        print(f"el error llega a 0.5 m a los {cifras(tiempo_para_error(b, 0.5), 4)} s  (4 cifras)")
    try:
        fig, ax = plt.subplots(figsize=(5.5, 3.2))
        ax.plot(ts, s_ver, color="black", ls="-"); ax.plot(ts, s_med, color="black", ls="--")
        rotula(ax, ts[-1], s_ver[-1], "verdadera", va="top", xytext=(4, -4))
        rotula(ax, ts[-1], s_med[-1], "medida (con sesgo)", va="bottom", xytext=(4, 4))
        ax.set_xlabel("t (s)"); ax.set_ylabel("s (m)")
        ax.set_title("posición verdadera (continua) y medida con el sesgo (discontinua)", fontsize=9); plt.show()
    except ERRORES:
        print("(no pude dibujar la gráfica con estos datos)")


widgets.interact(calculadora_movimiento,
    a_txt=widgets.Text(value="0.5 - 0.5*escalon(t - 2) - 0.5*escalon(t - 18) + 0.5*escalon(t - 20)",
                       description="a(t) =", continuous_update=False, layout=widgets.Layout(width="95%")),
    v0=widgets.FloatText(value=0, description="v0 (m/s)"),
    s0=widgets.FloatText(value=0, description="s0 (m)"),
    t_fin=widgets.FloatText(value=20, description="t final (s)"),
    b=widgets.FloatSlider(value=0.05, min=0, max=0.1, step=0.005, readout_format=".3f",
                          description="sesgo b", continuous_update=False));

# %%
# Casos de prueba de la sección 6.6 (resultado conocido)
_perfil = sp.Piecewise((sp.Rational(1, 2), x < 2), (0, x < 18), (-sp.Rational(1, 2), x < 20), (0, True))
_ts60 = np.linspace(0, 60, 601)

PRUEBAS_6 = [
    ("a = 2 desde el reposo: s(3) = 9 m", lambda: movimiento(2, 0, 0)[1].subs(x, 3) == 9),
    ("a = 6t desde el reposo: v = 3t², s(2) = 8 m",
     lambda: (lambda vs: sp.expand(vs[0] - 3 * x**2) == 0 and vs[1].subs(x, 2) == 8)(movimiento(6 * x, 0, 0))),
    ("perfil por tramos (Piecewise): s(20) = 18 m y v(20) = 0",
     lambda: (lambda vs: vs[1].subs(x, 20) == 18 and vs[0].subs(x, 20) == 0)(movimiento(_perfil, 0, 0))),
    ("el mismo perfil escrito con escalones: s(20) = 18 m",
     lambda: movimiento(parsear("0.5 - 0.5*escalon(t - 2) - 0.5*escalon(t - 18) + 0.5*escalon(t - 20)"))[1].subs(x, 20) == 18),
    ("trapecios con a = 2 en [0, 3]: s(3) = 9 m (exacto para a constante)",
     lambda: cerca(integra_datos(np.linspace(0, 3, 31), np.full(31, 2.0))[1][-1], 9, 1e-9)),
    ("sesgo 0.05 m/s² durante 60 s: error de 90 m (fórmula y trapecios)",
     lambda: cerca(error_sesgo(0.05, 60), 90) and cerca(integra_datos(_ts60, np.full(601, 0.05))[1][-1], 90, 1e-9)),
    ("con b = 0.05 el error llega a 0.5 m en √20 ≈ 4.472 s", lambda: cifras(tiempo_para_error(0.05, 0.5), 4) == "4.472"),
    ("tiempos que no aumentan se rechazan", lambda: _rechaza(lambda: integra_datos([0, 2, 1], [1, 1, 1]))),
]
for nombre, prueba in PRUEBAS_6:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 6.6).** Una piedra cae desde el reposo con $a=g=9.81$ m/s² (hacia abajo positivo). Calcula a mano cuánto cae en 2 s y escribe el número en metros.

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
        ref = movimiento(sp.Rational("9.81"), 0, 0)[1].subs(x, 2)
        print("La calculadora da: s(2) =", muestra(ref), "m")
        print("coinciden a 3 cifras" if coincide(mio, ref) else
              "NO coinciden: integra dos veces, v = g t y s = g t²/2; evalúa en t = 2 s")
