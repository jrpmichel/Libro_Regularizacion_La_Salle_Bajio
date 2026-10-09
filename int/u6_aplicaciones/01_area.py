# ID: INT-U6-NB01
# Notebook: int/u6_aplicaciones.ipynb · sección 6.1 área entre curvas
# Repositorio: int/u6_aplicaciones/01_area.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6.1 Área entre curvas
#
# $$A=\int_a^b\big|f(x)-g(x)\big|\,dx=\sum_k\left|\int_{c_k}^{c_{k+1}}\big(f(x)-g(x)\big)\,dx\right|,$$
#
# con $c_0=a$, $c_m=b$ y, entre ellos, los cruces de las curvas. Entre dos cruces $f-g$ no cambia de signo, así que cada parte se integra sin valor absoluto y se suma su magnitud. La integral con signo $\int_a^b(f-g)\,dx$ resta las partes donde $g$ queda arriba.
#
# En el ejemplo exacto de abajo, $y=4x$ y $y=x^3$ se cruzan en $x=2$ dentro de $[0,3]$.

# %%
h = 4 * x - x**3
print("cruces en (0, 3):", sp.solveset(h, x, sp.Interval.open(0, 3)))
I1, I2 = sp.integrate(h, (x, 0, 2)), sp.integrate(h, (x, 2, 3))
print("∫_0^2 (4x - x³) dx =", I1, "   ∫_2^3 (4x - x³) dx =", I2)
print("integral con signo =", I1 + I2, "   área =", abs(I1) + abs(I2))

# %% [markdown]
# **Calculadora de área entre curvas.** Escribe $f$, $g$ y los extremos ($a$ y $b$ aceptan `pi`, `pi/2` o fracciones). La calculadora busca los cruces dentro de $(a,b)$, integra $f-g$ en cada parte con su signo y suma las magnitudes. Si `sympy` no encuentra los cruces o la integral exacta, los calcula numéricamente y te lo dice. En la gráfica, el gris marca dónde $f$ queda arriba y el achurado dónde $g$ queda arriba.

# %%
def area_entre(f, g, a, b):
    """Área entre y = f(x) y y = g(x) en [a, b]. Devuelve un dict con los cruces dentro de (a, b), las partes
    [(p, q, integral con signo de f - g)], el área (suma de magnitudes), la integral con signo total y el método
    de integración. Un punto donde f - g solo toca 0 sin cambiar de signo no cuenta como cruce."""
    f, g = sp.sympify(f), sp.sympify(g)
    a, b = exacto(a, "a"), exacto(b, "b")
    if not a < b:
        raise ValueError("a debe ser menor que b.")
    h = f - g
    revisa_intervalo(h, a, b, "f - g")
    if simplifica(h) == 0:
        raise ValueError("f y g son la misma función: no encierran ninguna región (el área vale 0).")
    puntos, como = cruces(h, a, b)
    integ = Integrador()
    integ.solo_numerica = len(puntos) > 20          # f - g oscila mucho: todas las partes con quad
    partes = []
    for p, q in zip([a, *puntos], [*puntos, b]):
        I = integ(h, p, q)
        if partes and (float(partes[-1][2]) >= 0) == (float(I) >= 0):   # h solo tocó 0: se une con la parte anterior
            p0, _, I0 = partes.pop()
            p, I = p0, simplifica(I0 + I)
        partes.append((p, q, I))
    puntos = [q for _, q, _ in partes[:-1]]
    return {"cruces": puntos, "partes": partes, "area": simplifica(sum(abs(I) for _, _, I in partes)),
            "integral": simplifica(sum(I for _, _, I in partes)), "metodo": integ.metodo,
            "cruces_numericos": como == "numéricos" and bool(puntos)}


def calculadora_area(f_txt, g_txt, a_txt, b_txt):
    try:
        f, g = parsear(f_txt), parsear(g_txt)
        a, b = exacto(a_txt, "a"), exacto(b_txt, "b")
        r = area_entre(f, g, a, b)
    except ERRORES as err:
        print("Revisa la entrada:", explica(err)); return
    print(f"({CIFRAS} cifras significativas, redondeadas; si x y y están en metros, las áreas salen en m²)")
    print("cruces dentro de (a, b):", ", ".join(punto(c) for c in r["cruces"]) or "ninguno")
    if len(r["partes"]) > 12:
        print(f"   ({len(r['partes'])} partes; se muestran las primeras 12)")
    for p, q, I in r["partes"][:12]:
        arriba = "f arriba" if float(I) > 0 else "g arriba"
        print(f"   ∫ (f - g) dx de {punto(p)} a {punto(q)} = {muestra(I)} u²   ({arriba})")
    print(f"área total = {muestra(r['area'])} u²")
    print(f"integral con signo = {muestra(r['integral'])} u²   ({COMO[r['metodo']]}"
          + ("; cruces con brentq, numéricos" if r["cruces_numericos"] else "") + ")")
    try:
        xs = np.linspace(float(a), float(b), 600)
        yf, yg = para_graficar(f, xs), para_graficar(g, xs)
        fig, ax = plt.subplots(figsize=(5.5, 3.2))
        ax.fill_between(xs, yg, yf, where=yf >= yg, interpolate=True, color="0.8")
        ax.fill_between(xs, yg, yf, where=yf < yg, interpolate=True, facecolor="white", edgecolor="black", hatch="///")
        ax.plot(xs, yf, color="black", ls="-"); ax.plot(xs, yg, color="black", ls="--")
        rotula(ax, xs[-2], yf[-2], "f"); rotula(ax, xs[-2], yg[-2], "g")
        for c in r["cruces"][:50]:
            ax.axvline(float(c), color="0.4", ls=":", lw=0.8)
        ax.axhline(0, color="black", lw=0.5); ax.set_xlabel("x")
        ax.set_title("gris: f arriba de g · achurado: g arriba de f", fontsize=9); plt.show()
    except ERRORES:
        print("(no pude dibujar la gráfica con estos datos)")


widgets.interact(calculadora_area,
    f_txt=widgets.Text(value="4*x", description="f(x) =", continuous_update=False),
    g_txt=widgets.Text(value="x^3", description="g(x) =", continuous_update=False),
    a_txt=widgets.Text(value="0", description="a =", continuous_update=False),
    b_txt=widgets.Text(value="3", description="b =", continuous_update=False));

# %%
# Casos de prueba de la sección 6.1 (resultado conocido)
def _area_por_cuadratura(f, g, a, b):
    fn = numerica(f - g)
    return quad(lambda t: abs(float(fn(t))), float(a), float(b), limit=200)[0]


PRUEBAS_1 = [
    ("x + 2 y x² en [-1, 2]: área 9/2, sin cruces dentro",
     lambda: (lambda r: r["area"] == sp.Rational(9, 2) and r["cruces"] == [])(area_entre(x + 2, x**2, -1, 2))),
    ("4x y x³ en [0, 3]: cruce en 2, área 41/4, integral con signo -9/4",
     lambda: (lambda r: r["cruces"] == [2] and r["area"] == sp.Rational(41, 4) and r["integral"] == sp.Rational(-9, 4))(
         area_entre(4 * x, x**3, 0, 3))),
    ("sen x y cos x en [0, π]: cruce en π/4, área 2√2",
     lambda: (lambda r: r["cruces"] == [sp.pi / 4] and sp.simplify(r["area"] - 2 * sp.sqrt(2)) == 0)(
         area_entre(sp.sin(x), sp.cos(x), 0, sp.pi))),
    ("e^x y 3x en [0, 2]: dos cruces hallados numéricamente; el área coincide con ∫|f - g| por cuadratura",
     lambda: (lambda r: len(r["cruces"]) == 2 and r["cruces_numericos"]
              and cerca(r["area"], _area_por_cuadratura(sp.exp(x), 3 * x, 0, 2), 1e-8))(
         area_entre(sp.exp(x), 3 * x, 0, 2))),
    ("a ≥ b se rechaza", lambda: _rechaza(lambda: area_entre(x, x**2, 2, 1))),
    ("1/x y 0 en [-1, 1] se rechaza (asíntota en 0)", lambda: _rechaza(lambda: area_entre(1 / x, 0, -1, 1))),
]
for nombre, prueba in PRUEBAS_1:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 6.1).** Calcula a mano el área de la región encerrada entre $y=x^2$ y $y=2x$ y escribe el número.

# %%
mi_valor = None     # escribe un número (también sirve "4/3")

if mi_valor is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    try:
        mio = a_numero(mi_valor)
    except ValueError as err:
        print("Revisa tu número:", err)
    else:
        ref = area_entre(2 * x, x**2, 0, 2)["area"]
        print("La calculadora da:", muestra(ref), "u²")
        print("coinciden a 3 cifras" if coincide(mio, ref) else
              "NO coinciden: las curvas se cortan en x = 0 y x = 2; integra 2x - x² entre esos valores")
