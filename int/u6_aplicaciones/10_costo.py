# ID: INT-U6-NB10
# Notebook: int/u6_aplicaciones.ipynb · sección 6.10 costo acumulado
# Repositorio: int/u6_aplicaciones/10_costo.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6.10 Costo acumulado
#
# $$C(q_1)-C(q_0)=\int_{q_0}^{q_1}C'(q)\,dq,\qquad C(T)=P+\int_0^Tm(t)\,dt.$$
#
# El costo marginal $C'(q)$ es lo que cuesta producir la pieza número $q$; integrarlo da el costo de un lote. Para una máquina con precio de compra $P$ y tasa de mantenimiento $m(t)$ en pesos por año, $C(T)$ es lo que llevas gastado a los $T$ años. Entre dos máquinas, la más cara al comprar puede salir más barata con el tiempo si su mantenimiento crece más despacio; el punto de equilibrio es el $T>0$ donde $C_A(T)=C_B(T)$.
#
# En el ejemplo exacto, la máquina A cuesta 250 000 pesos con $m_A=20\,000+6\,000t$ y la B cuesta 400 000 pesos con $m_B=12\,000+1\,500t$.

# %%
T_s, t_m = sp.symbols("T t", nonnegative=True)
C_A = 250000 + sp.integrate(20000 + 6000 * t_m, (t_m, 0, T_s))
C_B = 400000 + sp.integrate(12000 + 1500 * t_m, (t_m, 0, T_s))
print("C_A(T) =", C_A, "   C_B(T) =", C_B)
print("equilibrio:", sp.solveset(C_A - C_B, T_s, sp.Interval.open(0, sp.oo)), "años")

# %% [markdown]
# **Calculadora de equilibrio entre dos máquinas.** Escribe el precio de compra de cada máquina en pesos y su tasa de mantenimiento $m(t)$ en pesos por año, con $t$ en años. La calculadora busca los tiempos positivos en que los costos acumulados se igualan (exactos con `sympy` si se puede; si no, numéricos dentro del horizonte) y dibuja $C_A(T)$ y $C_B(T)$.

# %%
_TAU10 = sp.Dummy("tau", real=True)


def costo_total(C_marginal, q0, q1, fijo=0):
    """Costo de producir de q0 a q1: fijo + ∫_{q0}^{q1} C'(q) dq (en las unidades de C' por unidad de q)."""
    Cm = sp.sympify(C_marginal)
    q0, q1, fijo = exacto(q0, "q0"), exacto(q1, "q1"), exacto(fijo, "el costo fijo")
    if q0 < 0 or not q0 < q1:
        raise ValueError("necesitas 0 ≤ q0 < q1.")
    revisa_intervalo(Cm, q0, q1, "C'", var="q")
    return simplifica(fijo + integra(Cm, q0, q1)[0])


def costo_acumulado(precio, m, T):
    """C(T) = precio + ∫_0^T m(t) dt: lo gastado en una máquina a los T años (m en pesos por año)."""
    precio, T = exacto(precio, "el precio"), exacto(T, "T")
    if T < 0:
        raise ValueError("T es un tiempo: no puede ser negativo.")
    if T == 0:
        return precio
    m = sp.sympify(m)
    revisa_intervalo(m, 0, T, "m", var="t")
    return simplifica(precio + integra(m, 0, T)[0])


def _acumulada(m):
    """∫_0^x m(t) dt como expresión en x (para resolver C_A = C_B), o None si sympy no da una fórmula útil."""
    m = sp.sympify(m)
    if m.has(sp.Heaviside):
        m = m.rewrite(sp.Piecewise)
    try:
        r = con_limite(lambda: sp.integrate(m.subs(x, _TAU10), (_TAU10, 0, x)), 10)
    except (ValueError, NotImplementedError, TypeError, ArithmeticError):
        return None
    if r.has(*_RARAS) or r.has(sp.I, sp.oo, -sp.oo, sp.zoo, sp.nan):
        return None
    return r


def _acumula_malla(f, ts):
    """∫_0^t f en cada punto de la malla, con la regla del punto medio (no evalúa f en t = 0, donde puede no existir)."""
    v = f((ts[1:] + ts[:-1]) / 2)
    return np.concatenate(([0.0], np.cumsum(np.diff(ts) * v)))


def equilibrio(precio_A, m_A, precio_B, m_B, horizonte=100, con_metodo=False):
    """Tiempos T > 0 (años) en que C_A(T) = C_B(T), en orden. Exactos con sympy si se puede (todos los T > 0);
    si no, se buscan numéricamente solo en (0, horizonte]. Lista vacía si no se igualan. Con con_metodo=True
    devuelve (tiempos, "exacta" o "numérica")."""
    pA, pB = exacto(precio_A, "el precio de A"), exacto(precio_B, "el precio de B")
    m_A, m_B = sp.sympify(m_A), sp.sympify(m_B)
    if not horizonte > 0:
        raise ValueError("el horizonte debe ser positivo.")
    for m, nombre in ((m_A, "m_A"), (m_B, "m_B")):
        revisa_intervalo(m, 0, horizonte, nombre, var="t")
    def _sale(tiempos, metodo):
        return (tiempos, metodo) if con_metodo else tiempos
    acum = _acumulada(m_A - m_B)
    if acum is not None:
        D = (pA - pB) + acum
        if simplifica(D) == 0:
            raise ValueError("las dos máquinas cuestan lo mismo en todo momento: no hay un punto de equilibrio único.")
        try:
            sol = con_limite(lambda: sp.solveset(D, x, sp.Interval.open(0, sp.oo)), 10)
        except (ValueError, NotImplementedError, TypeError):
            sol = None
        if sol == sp.EmptySet:
            return _sale([], "exacta")
        if isinstance(sol, sp.FiniteSet) and all(c.is_real for c in sol):
            return _sale(sorted((simplifica(c) for c in sol), key=float), "exacta")
    fd = numerica(m_A - m_B)                        # numérico: cambios de signo de C_A - C_B en una malla
    ts = np.linspace(0, float(horizonte), 4001)
    D = float(pA - pB) + _acumula_malla(fd, ts)
    def D_exacta(t):
        return float(pA - pB) + quad(lambda s_: float(fd(s_)), 0, t, limit=200)[0]
    raices = []
    for k in range(1, len(ts) - 1):
        if D[k] == 0:
            raices.append(ts[k])
        elif D[k] * D[k + 1] < 0:
            try:
                raices.append(brentq(D_exacta, ts[k], ts[k + 1], xtol=1e-12))
            except ValueError:                      # quad y la malla difieren en el signo: interpolación lineal
                raices.append(ts[k] - D[k] * (ts[k + 1] - ts[k]) / (D[k + 1] - D[k]))
    return _sale([sp.Float(r, 12) for r in raices], "numérica")


def calculadora_equilibrio(precio_A, m_A_txt, precio_B, m_B_txt, horizonte):
    try:
        m_A, m_B = parsear(m_A_txt), parsear(m_B_txt)
        if precio_A < 0 or precio_B < 0:
            raise ValueError("los precios no pueden ser negativos.")
        tiempos, metodo = equilibrio(precio_A, m_A, precio_B, m_B, horizonte, con_metodo=True)
        costos = [costo_acumulado(precio_A, m_A, T) for T in tiempos]
    except ERRORES as err:
        print("Revisa la entrada:", explica(err)); return
    exacta = metodo == "exacta"
    print(f"({CIFRAS} cifras significativas, redondeadas; T en años, costos en pesos; "
          + ("cálculo exacto con sympy, para todo T > 0)" if exacta else
             f"búsqueda numérica, solo en (0, {horizonte}] años)"))
    if not tiempos:
        print("los costos acumulados no se igualan para ningún T > 0" if exacta else
              f"los costos acumulados no se igualan en (0, {horizonte}] años; prueba un horizonte mayor")
    for T, C in zip(tiempos, costos):
        print(f"equilibrio en T = {muestra(T)} años: los dos costos llegan a {cifras(C, CIFRAS)} pesos"
              + ("   (fuera del horizonte de la gráfica)" if float(T) > horizonte else ""))
    try:
        def _barata(T):
            dif = float(costo_acumulado(precio_A, m_A, T) - costo_acumulado(precio_B, m_B, T))
            return "A" if dif < 0 else ("B" if dif > 0 else "ninguna")
        if tiempos:
            print(f"antes del primer equilibrio sale más barata la máquina {_barata(float(tiempos[0]) / 2)}; "
                  f"después del último, la {_barata(float(tiempos[-1]) + 1)}")
        Ts = np.linspace(0, horizonte, 600)
        CA = precio_A + _acumula_malla(numerica(m_A), Ts)
        CB = precio_B + _acumula_malla(numerica(m_B), Ts)
        fig, ax = plt.subplots(figsize=(5.5, 3.2))
        ax.plot(Ts, CA, color="black", ls="-"); ax.plot(Ts, CB, color="black", ls="--")
        rotula(ax, Ts[-1], CA[-1], "A"); rotula(ax, Ts[-1], CB[-1], "B")
        for T, C in zip(tiempos, costos):
            if float(T) <= horizonte:
                ax.axvline(float(T), color="0.4", ls=":", lw=1)
                rotula(ax, float(T), float(C), f"T = {cifras(T, CIFRAS)} años", xytext=(-6, 12), ha="right")
        ax.set_xlabel("T (años)"); ax.set_ylabel("costo acumulado (pesos)")
        ax.ticklabel_format(axis="y", style="plain"); plt.show()
    except ERRORES:
        print("(no pude dibujar la gráfica con estos datos)")


_ancho = {"description_width": "initial"}
widgets.interact(calculadora_equilibrio,
    precio_A=widgets.FloatText(value=250000, description="precio A (pesos)", style=_ancho),
    m_A_txt=widgets.Text(value="20000 + 6000*t", description="m_A(t) (pesos/año) =", continuous_update=False, style=_ancho),
    precio_B=widgets.FloatText(value=400000, description="precio B (pesos)", style=_ancho),
    m_B_txt=widgets.Text(value="12000 + 1500*t", description="m_B(t) (pesos/año) =", continuous_update=False, style=_ancho),
    horizonte=widgets.IntSlider(value=12, min=1, max=30, description="horizonte (años)", style=_ancho,
                                continuous_update=False));

# %%
# Casos de prueba de la sección 6.10 (resultado conocido)
PRUEBAS_10 = [
    ("C'(q) = 120 - 0.4q de q = 100 a 150: 3500 pesos", lambda: costo_total(120 - sp.Rational(2, 5) * x, 100, 150) == 3500),
    ("máquinas A y B: equilibrio en T = (-16 + 2√1414)/9 ≈ 6.5785 años",
     lambda: (lambda c: len(c) == 1 and sp.simplify(c[0] - (-16 + 2 * sp.sqrt(1414)) / 9) == 0 and cifras(c[0], 5) == "6.5785")(
         equilibrio(250000, 20000 + 6000 * x, 400000, 12000 + 1500 * x))),
    ("C_A(5) = 425000 pesos", lambda: costo_acumulado(250000, 20000 + 6000 * x, 5) == 425000),
    ("C_B(10) = 595000 pesos", lambda: costo_acumulado(400000, 12000 + 1500 * x, 10) == 595000),
    ("curva de aprendizaje 10·q^(-0.32) de 0 a 100 piezas: 336.89 minutos",
     lambda: round(float(costo_total(parsear("10*q^(-0.32)"), 0, 100)), 2) == 336.89),
    ("mismo mantenimiento y distinto precio: nunca se igualan", lambda: equilibrio(100, 1, 200, 1) == []),
    ("q1 < q0 se rechaza", lambda: _rechaza(lambda: costo_total(10, 50, 0))),
]
for nombre, prueba in PRUEBAS_10:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 6.10).** Cada pieza cuesta 10 pesos, sin importar cuántas hagas ($C'(q)=10$). Calcula a mano el costo de producir las primeras 50 piezas, sin costo fijo, y escribe el número en pesos.

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
        ref = costo_total(10, 0, 50)
        print("La calculadora da:", muestra(ref), "pesos")
        print("coinciden a 3 cifras" if coincide(mio, ref) else
              "NO coinciden: ∫_0^50 10 dq es el área de un rectángulo de base 50 y altura 10")
