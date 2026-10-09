# ID: INT-U5-NB05
# Notebook: int/u5_tecnicas.ipynb · sección 5.5 integrales impropias
# Repositorio: int/u5_tecnicas/05_impropias.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 5. Integrales impropias
#
# $$\int_a^{\infty}f(x)\,dx=\lim_{R\to\infty}\int_a^Rf(x)\,dx,\qquad\int_a^bf(x)\,dx=\lim_{\varepsilon\to0^+}\int_{a+\varepsilon}^bf(x)\,dx\ \ (f\text{ no acotada cerca de }a).$$
#
# Si hay un punto problemático dentro del intervalo, o los dos extremos son infinitos, la integral se separa, y **cada parte debe converger por su lado**. Integrales $p$: $\int_1^{\infty}x^{-p}dx$ converge si $p>1$; $\int_0^1x^{-p}dx$ converge si $p<1$.

# %%
print(" R       ∫_1^R dx/x²   ∫_1^R dx/x")
for R in (10, 100, 1000, 10000):
    print(f"{R:<7d} {cifras(1 - 1 / R, 6):>12} {cifras(math.log(R), 6):>12}")

# %% [markdown]
# **Explorador de integrales impropias.** Escribe $f(x)$ y los extremos (pueden ser `-inf`, `inf`, `pi/2` o una fracción). El explorador busca los puntos de $[a,b]$ donde $f$ no está definida, separa la integral, arma la tabla de integrales parciales de cada parte y da un veredicto por parte y para la integral completa. El veredicto sale de un cálculo exacto (`sympy.integrate` o el límite de la antiderivada); si solo hay evidencia numérica, el explorador lo dice y no afirma que la integral converja o diverja. Si una parte diverge, no suma las demás.

# %%
def lee_extremo(texto):
    """'inf', '-inf', 'pi/2', '3', '0.5' -> número de sympy (exacto si se puede)."""
    t = str(texto).strip().lower().replace(" ", "").replace("∞", "inf")
    if t in ("inf", "+inf", "oo", "+oo", "infinito"):
        return sp.oo
    if t in ("-inf", "-oo", "-infinito"):
        return -sp.oo
    try:
        v = parsear(t)
    except ValueError as err:
        raise ValueError("Escribe los extremos como números (3, 0.5, 1/2, pi/2), inf o -inf.") from err
    if v.free_symbols:
        raise ValueError("Los extremos son números: no lleves x en ellos.")
    return v


def _txt(v):
    """Extremo legible: ∞, un entero o fracción tal cual, o el valor exacto con su aproximación."""
    if v in (sp.oo, -sp.oo):
        return "∞" if v == sp.oo else "-∞"
    return str(v) if v.is_Rational else f"{v} ≈ {cifras(v, 6)}"


def _problemas(f, a, b):
    """Puntos de [a, b] donde f no es continua; ValueError si no son un número finito de puntos."""
    try:
        dom = con_limite(lambda: sp.calculus.util.continuous_domain(f, x, sp.Interval(a, b)), 10)
    except (ValueError, NotImplementedError):
        raise ValueError("no pude determinar dónde f no está definida; prueba con un intervalo más corto.") from None
    fuera = sp.Interval(a, b) - dom
    if fuera == sp.EmptySet:
        return []
    if isinstance(fuera, sp.FiniteSet):
        return sorted(fuera, key=float)
    if isinstance(fuera, sp.ConditionSet) or fuera.has(sp.ConditionSet):
        raise ValueError("no pude determinar dónde f no está definida; prueba con un intervalo más corto.")
    if fuera.is_iterable or fuera.has(sp.ImageSet):
        raise ValueError("f tiene una infinidad de puntos problemáticos en el intervalo (como tan x en [0, ∞)): "
                         "el explorador no puede separar la integral en un número finito de partes.")
    raise ValueError("f no está definida (o no es real) en todo un tramo del intervalo.")


def _revisa_real(f, p, q):
    """ValueError si f da valores complejos dentro de la parte (p, q) (por ejemplo x^(-1/3) con x < 0)."""
    lo = float(p) if p != -sp.oo else min(float(q) - 10, -10.0)
    hi = float(q) if q != sp.oo else max(float(p) + 10, 10.0)
    for k in (0.17, 0.41, 0.63, 0.88):
        try:
            valor_real(f, sp.nsimplify(lo + k * (hi - lo)))
        except ValueError:
            raise ValueError(f"f no es real dentro de [{_txt(p)}, {_txt(q)}]. Si escribiste una raíz cúbica de un "
                             "número negativo, sympy la toma compleja; usa un intervalo donde la base sea positiva.") from None


def _partes(a, b, malos):
    """Separa [a, b] para que cada parte tenga a lo más un extremo problemático."""
    cortes = sorted({a, b, *malos}, key=lambda v: float(v))
    partes = []
    for p, q in zip(cortes, cortes[1:]):
        mp, mq = (p == -sp.oo or p in malos), (q == sp.oo or q in malos)
        if mp and mq:
            if p == -sp.oo and q == sp.oo:
                c = sp.Integer(0)
            elif q == sp.oo:
                c = p + 1
            elif p == -sp.oo:
                c = q - 1
            else:
                c = (p + q) / 2
            partes += [(p, c, "izq"), (c, q, "der")]
        else:
            partes.append((p, q, "izq" if mp else ("der" if mq else None)))
    return partes


def _quad(fn, u, v):
    """∫_u^v por cuadratura, en tramos de longitud <= 50 (las oscilaciones largas no confunden a quad)."""
    cortes = list(np.linspace(u, v, max(2, int(math.ceil(abs(v - u) / 50)) + 1)))
    total = 0.0
    with np.errstate(all="ignore"):
        for s_, t_ in zip(cortes, cortes[1:]):
            total += quad(fn, s_, t_, limit=200)[0]
    return total


def _acercamiento(p, q, lado):
    """Los cuatro puntos z de la tabla, cada vez más cerca del extremo problemático."""
    largo = 1.0 if (p == -sp.oo or q == sp.oo) else float(q - p)
    eps = [10.0**-k if largo >= 1 else largo / 10**k for k in (1, 2, 3, 4)]
    if lado == "der":
        return [(10**k if 10**k > p else p + 10**k) for k in (1, 2, 3, 4)] if q == sp.oo else [float(q) - e for e in eps]
    return [(-10**k if -10**k < q else q - 10**k) for k in (1, 2, 3, 4)] if p == -sp.oo else [float(p) + e for e in eps]


def _parciales(F, f, p, q, lado):
    """Tabla de integrales parciales hacia el extremo problemático. Con antiderivada elemental, F(z) - F(p);
    sin ella, sumas de cuadraturas que se acumulan tramo por tramo desde el extremo bueno."""
    zs, filas = _acercamiento(p, q, lado), []
    if F is not None:
        for z in zs:
            try:
                v = valor_real(F, z) - valor_real(F, p) if lado == "der" else valor_real(F, q) - valor_real(F, z)
            except OverflowError:
                v = math.inf
            except (ValueError, TypeError, ZeroDivisionError):
                v = float("nan")
            filas.append((z, v))
        return filas
    fn = sp.lambdify(x, f, "math")
    bueno = float(p) if lado == "der" else float(q)
    acumulado, previo = 0.0, bueno
    for z in zs:
        try:
            tramo = _quad(fn, previo, float(z))
            acumulado += tramo if lado == "der" else -tramo
        except OverflowError:
            acumulado = math.inf
        except (ValueError, ZeroDivisionError):
            acumulado = float("nan")
        previo = float(z)
        filas.append((z, acumulado))
    return filas


def _veredicto(I):
    """(valor, estado) a partir de un resultado exacto de sympy, o None si no sirve para decidir."""
    if I is None or I.has(sp.Integral, sp.nan) or I is sp.nan:
        return None
    if isinstance(I, sp.AccumBounds) or I.has(sp.AccumBounds):
        return None, "oscila"
    if I.has(sp.oo, -sp.oo, sp.zoo) or I.is_finite is False:
        return None, "diverge"
    if I.is_real is False or I.has(sp.I):
        return None
    return sp.simplify(I), "converge"


def _valor_parte(f, F, p, q, lado):
    """(valor, estado) de una parte: estado es converge, diverge, oscila o no concluyente."""
    try:
        r = _veredicto(con_limite(lambda: sp.integrate(f, (x, p, q)), 15))
    except (ValueError, NotImplementedError):
        r = None
    if r is None and F is not None and lado is not None:
        try:
            if lado == "der":
                lim = con_limite(lambda: sp.limit(F, x, q, "-"), 10)
                r = _veredicto(lim - F.subs(x, p))
            else:
                lim = con_limite(lambda: sp.limit(F, x, p, "+"), 10)
                r = _veredicto(F.subs(x, q) - lim)
        except (ValueError, NotImplementedError):
            r = None
    return r if r is not None else (None, "no concluyente")


def explorar(f, a, b):
    """Lista de partes (p, q, tabla, valor, estado) y el veredicto global (valor, estado)."""
    if not (a < b):
        raise ValueError("el extremo izquierdo debe ser menor que el derecho.")
    malos = _problemas(f, a, b)
    try:
        F = antiderivada(f)
    except ValueError:
        F = None
    resultado = []
    for p, q, lado in _partes(a, b, malos):
        _revisa_real(f, p, q)
        valor, estado = _valor_parte(f, F, p, q, lado)
        tabla = _parciales(F, f, p, q, lado) if lado is not None else []
        resultado.append((p, q, tabla, valor, estado))
    estados = [r[4] for r in resultado]
    if "diverge" in estados or "oscila" in estados:
        return resultado, (None, "diverge")
    if "no concluyente" in estados:
        return resultado, (None, "no concluyente")
    return resultado, (sp.simplify(sum(r[3] for r in resultado)), "converge")


_LEYENDA = {"diverge": "DIVERGE", "oscila": "DIVERGE: las integrales parciales oscilan sin acercarse a un número",
            "no concluyente": "NO CONCLUYENTE: sympy no dio un valor exacto; mira la tabla, pero una tabla sola "
                              "no demuestra que la integral converja o diverja"}


def explorador_impropias(f_txt, a_txt, b_txt):
    try:
        f = parsear(f_txt)
        a, b = lee_extremo(a_txt), lee_extremo(b_txt)
        partes, (total, estado) = explorar(f, a, b)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print("(ángulos en radianes; parciales con 6 cifras significativas, redondeadas)")
    for p, q, tabla, valor, est in partes:
        print(f"parte {'(' if p == -sp.oo else '['}{_txt(p)}, {_txt(q)}{')' if q == sp.oo else ']'}:")
        for z, v in tabla:
            print(f"   hasta {cifras(z, 6):>10}:  {'más de 1e308' if math.isinf(v) else cifras(v, 6)}")
        print("   →", f"converge a {valor} ≈ {cifras(valor, 6)}" if est == "converge" else _LEYENDA[est])
    if estado == "diverge":
        print("La integral completa DIVERGE: al menos una parte diverge, y las partes divergentes no se suman.")
    elif estado == "no concluyente":
        print("La integral completa queda SIN VEREDICTO: al menos una parte solo tiene evidencia numérica.")
    else:
        print(f"La integral completa converge a {total} ≈ {cifras(total, 6)}")


widgets.interact(explorador_impropias,
    f_txt=widgets.Text(value="1/x^2", description="f(x) =", continuous_update=False),
    a_txt=widgets.Text(value="1", description="a =", continuous_update=False),
    b_txt=widgets.Text(value="inf", description="b =", continuous_update=False));

# %%
# Casos de prueba de la sección 5 (los del prompt del libro y cuatro más)
def _valor(f, a, b):
    return explorar(f, a, b)[1][0]


def _estado(f, a, b):
    return explorar(f, a, b)[1][1]


PRUEBAS_5 = [
    ("e^(-2x) en [0, ∞): 1/2", lambda: _valor(sp.exp(-2 * x), 0, sp.oo) == sp.Rational(1, 2)),
    ("1/x en [1, ∞): diverge; parciales 2.30259, 4.60517, 6.90776, 9.21034",
     lambda: (lambda r: r[1][1] == "diverge" and [cifras(v, 6) for _, v in r[0][0][2]] == ["2.30259", "4.60517", "6.90776", "9.21034"])(
         explorar(1 / x, 1, sp.oo))),
    ("1/√x en [0, 1]: 2", lambda: _valor(1 / sp.sqrt(x), 0, 1) == 2),
    ("1/x en [-1, 1]: diverge aunque sea impar", lambda: _estado(1 / x, -1, 1) == "diverge"),
    ("1/(1 + x²) en (-∞, ∞): π", lambda: _valor(1 / (1 + x**2), -sp.oo, sp.oo) == sp.pi),
    ("x^(-1/3) en [0, 1]: 3/2", lambda: _valor(x**sp.Rational(-1, 3), 0, 1) == sp.Rational(3, 2)),
    ("sen x / x en [0, ∞): π/2 (sin antiderivada elemental)", lambda: sp.simplify(_valor(sp.sin(x) / x, 0, sp.oo) - sp.pi / 2) == 0),
    ("sen x en [0, ∞): diverge porque oscila", lambda: (lambda r: r[0][0][4] == "oscila" and r[1][1] == "diverge")(
        explorar(sp.sin(x), 0, sp.oo))),
]
for nombre, prueba in PRUEBAS_5:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 5).** Calcula a mano $\int_1^{\infty}\frac{dx}{x^3}$ y escribe el número.

# %%
mi_valor = None     # escribe un número (también sirve "1/2")

if mi_valor is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    try:
        mio = a_numero(mi_valor)
    except ValueError as err:
        print("Revisa tu número:", err)
    else:
        ref = _valor(1 / x**3, 1, sp.oo)
        print("El explorador da:", ref)
        print("coinciden a 3 cifras" if coincide(mio, ref) else
              "NO coinciden: ∫_1^R x^(-3) dx = 1/2 - 1/(2R²); toma el límite")
