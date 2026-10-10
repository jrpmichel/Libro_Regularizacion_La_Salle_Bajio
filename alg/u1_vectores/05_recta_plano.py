# ID: ALG-U1-NB05
# Notebook: alg/u1_vectores.ipynb · sección 1.5 rectas y planos en el espacio
# Repositorio: alg/u1_vectores/05_recta_plano.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 5. Rectas y planos en el espacio
#
# * **Recta** por $P_0(x_0,y_0,z_0)$ con dirección $\mathbf u=(a,b,c)\neq\mathbf 0$: ecuación vectorial $\mathbf r(t)=P_0+t\,\mathbf u$ y paramétricas $x=x_0+at,\;y=y_0+bt,\;z=z_0+ct$, con $t$ real.
# * **Plano** por $P_0$ con normal $\mathbf n=(a,b,c)\neq\mathbf 0$: los puntos $P$ con $\mathbf n\cdot(P-P_0)=0$. Al desarrollar queda la ecuación cartesiana $ax+by+cz=d$, con $d=\mathbf n\cdot P_0$.
# * **Plano por tres puntos** $A$, $B$, $C$: una normal es $\mathbf n=\overrightarrow{AB}\times\overrightarrow{AC}$. Si sale $\mathbf 0$, los puntos están alineados y no definen un plano único.
# * **Distancia** de un punto $Q$ al plano: $D=\dfrac{|\mathbf n\cdot Q-d|}{\|\mathbf n\|}$. El signo de $\mathbf n\cdot Q-d$ dice de qué lado queda $Q$: positivo, del lado hacia donde apunta $\mathbf n$.
# * **Intersección** de la recta con el plano: al sustituir $\mathbf r(t)$ queda $\mathbf n\cdot P_0+t\,(\mathbf n\cdot\mathbf u)=d$.
#
# | Caso | Condición | Resultado |
# |---|---|---|
# | corta en un punto | $\mathbf n\cdot\mathbf u\neq0$ | $t^*=\dfrac{d-\mathbf n\cdot P_0}{\mathbf n\cdot\mathbf u}$, en el punto $P_0+t^*\mathbf u$ |
# | paralela | $\mathbf n\cdot\mathbf u=0$ y $\mathbf n\cdot P_0\neq d$ | no toca el plano |
# | contenida | $\mathbf n\cdot\mathbf u=0$ y $\mathbf n\cdot P_0=d$ | toda la recta está en el plano |
#
# Cualquier múltiplo no nulo de $\mathbf n$ (con $d$ multiplicado por lo mismo) describe el mismo plano. Con tres puntos, la calculadora simplifica: si las coordenadas son enteras o fracciones (como `7/5` o `0.25`), hace las cuentas con fracciones exactas y da coeficientes enteros sin factor común; si no (por ejemplo, con `sqrt(2)`), los divide entre el de mayor valor absoluto. Además deja positivo el primer coeficiente que no es cero. Si $t^*<0$, el corte queda antes de $P_0$: un rayo (solo $t\ge0$), como el de una cámara, no llega a él.

# %%
def _normal(n):
    n = vector(n, 3, "la normal n")
    if not np.any(n):
        raise ValueError("la normal n no puede ser el vector cero.")
    return n


def _enteros(valores):
    """True si todos son enteros (dentro del redondeo). 1e-200 no cuenta como el entero 0."""
    return all(abs(v) < 2**53 and abs(v - round(v)) <= 1e-9 * max(1.0, abs(v)) and (v == 0 or round(v) != 0)
               for v in valores)


def _num_txt(v, n=4):
    """Entero exacto sin decimales; lo demás, redondeado a n cifras significativas y sin ceros de relleno que
    parezcan exactos: 123456789.4 con 4 cifras se escribe 1.235×10^8, no 123500000."""
    if _enteros([v]):
        return str(int(round(v)))
    texto = cifras(v, n)
    if "." in texto or "×" in texto or len(texto.lstrip("-")) <= n:
        return texto
    q = Decimal(texto)                            # entero con ceros de relleno: a notación científica
    e = q.adjusted()
    return f"{format(q.scaleb(-e).quantize(Decimal(1).scaleb(1 - n)), 'f')}×10^{e}"


def _fmt_txt(v, n=4):
    """Vector '(a, b, c)' con el mismo formato que los coeficientes del plano."""
    return "(" + ", ".join(_num_txt(c, n) for c in v) + ")"


def _lado(P, n, d):
    """n · P − d, con la cancelación por redondeo limpiada."""
    terminos = [*(n * P), -d]
    s = math.fsum(terminos)
    return 0.0 if abs(s) <= 1e-14 * math.fsum(abs(x) for x in terminos) else s


def _fraccion(x):
    """x como fracción de denominador hasta 10^6 si x es exactamente el flotante de esa fracción (7/5, 0.25,
    -9/7); si no (sqrt(2), pi), None."""
    f = Fraction(x).limit_denominator(10**6)
    return f if float(f) == x else None


def _plano_exacto(A, B, C):
    """(n, d) con enteros sin factor común, calculados con fracciones exactas; None si alguna coordenada
    no es una fracción sencilla o los enteros no caben en un flotante."""
    fr = [_fraccion(x) for x in (*A, *B, *C)]
    if any(f is None for f in fr):
        return None
    a, b, c = fr[0:3], fr[3:6], fr[6:9]
    u = [q - p for p, q in zip(a, b)]
    w = [q - p for p, q in zip(a, c)]
    n = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]]
    if not any(n):
        return None
    d = sum(ni * ai for ni, ai in zip(n, a))
    m = math.lcm(*(f.denominator for f in (*n, d)))
    enteros = [int(f * m) for f in (*n, d)]
    g = math.gcd(*enteros)
    try:
        return np.array([float(e // g) for e in enteros[:3]]), float(enteros[3] // g)
    except OverflowError:
        return None


def plano_tres_puntos(A, B, C):
    """(n, d) del plano n · (x, y, z) = d que pasa por A, B y C, con n = AB × AC simplificado: enteros sin
    factor común si las coordenadas son enteras o fracciones sencillas; si no, divididos entre el coeficiente
    de mayor valor absoluto. El primer coeficiente distinto de cero queda positivo."""
    A, B, C = vector(A, 3, "A"), vector(B, 3, "B"), vector(C, 3, "C")
    AB, AC = B - A, C - A
    if paralelos(AB, AC):
        raise ValueError("los puntos están alineados (o repetidos): no definen un plano único.")
    exacto = _plano_exacto(A, B, C)
    if exacto is not None:
        n, d = exacto
    else:
        n, _ = cruz(AB, AC)
        d = _punto(n, A)
        escala = float(np.max(np.abs(n)))
        n, d = n / escala, d / escala
    if n[np.flatnonzero(n)[0]] < 0:
        n, d = -n, -d
    return n + 0.0, float(d) + 0.0


def plano_texto(n, d, n_cifras=4):
    """'3x + 5y + 4z = 13' a partir de la normal n y de d."""
    n, d = _normal(n), numero(d, "d")
    texto = ""
    for c, var in zip(n, "xyz"):
        if c == 0:
            continue
        coef = _num_txt(abs(c), n_cifras)
        termino = ("" if coef == "1" else coef) + var
        texto += (("-" if c < 0 else "") + termino) if not texto else ((" - " if c < 0 else " + ") + termino)
    return f"{texto} = {_num_txt(d, n_cifras)}"


def distancia(P, n, d):
    """(D, lado): distancia del punto P al plano n · (x, y, z) = d y de qué lado queda."""
    P, n, d = vector(P, 3, "el punto"), _normal(n), numero(d, "d")
    s = _lado(P, n, d)
    lado = ("sobre el plano" if s == 0 else
            "del lado hacia donde apunta n" if s > 0 else "del lado contrario a n")
    return abs(s) / _norma(n), lado


def interseccion(P0, u, n, d):
    """('punto', t*, X), ('paralela', None, None) o ('contenida', None, None) para la recta P0 + t u."""
    P0, u, n, d = vector(P0, 3, "P0"), vector(u, 3, "la dirección u"), _normal(n), numero(d, "d")
    if not np.any(u):
        raise ValueError("la dirección u de la recta no puede ser el vector cero.")
    s = _lado(P0, n, d)
    den = _punto(n, u)
    if abs(den) <= 1e-12 * _norma(n) * _norma(u):          # n · u = 0: la recta va paralela al plano
        return ("contenida", None, None) if s == 0 else ("paralela", None, None)
    t = -s / den + 0.0
    return "punto", t, _limpia(P0 + t * u, np.abs(P0) + np.abs(t * u))


def misma_ecuacion(p, q, tol=1e-9):
    """True si (a, b, c, d) y (a', b', c', d') son múltiplos no nulos uno del otro (el mismo plano)."""
    p, q = vector(p, 4, "la ecuación (a, b, c, d)"), vector(q, 4, "la ecuación (a, b, c, d)")
    if not np.any(p[:3]) or not np.any(q[:3]):
        raise ValueError("a, b y c no pueden ser los tres 0: la normal sería el vector cero.")
    up, uq = p / _norma(p), q / _norma(q)
    return min(_norma(up - uq), _norma(up + uq)) <= tol


def _parametricas(P0, u, n_cifras):
    """'x = 1 + 2t,   y = -t,   z = 3' (sin términos que valen 0)."""
    partes = []
    for x0, a, var in zip(P0, u, "xyz"):
        coef = _num_txt(abs(a), n_cifras)
        termino = ("" if coef == "1" else coef) + "t"
        if a == 0:
            partes.append(f"{var} = {_num_txt(x0, n_cifras)}")
        elif x0 == 0:
            partes.append(f"{var} = {'-' if a < 0 else ''}{termino}")
        else:
            partes.append(f"{var} = {_num_txt(x0, n_cifras)} {'-' if a < 0 else '+'} {termino}")
    return ",   ".join(partes)


def _pie(P, n, d):
    """Punto del plano más cercano a P (pie de la perpendicular), con la normal unitaria: sin desbordes."""
    nn = _norma(n)
    return P - (_lado(P, n, d) / nn) * (n / nn)


def _base_plano(n):
    """Dos vectores unitarios perpendiculares entre sí y a n (para dibujar el plano)."""
    nu = n / _norma(n)
    e1 = np.cross(nu, np.eye(3)[np.argmin(np.abs(nu))])
    e1 = e1 / _norma(e1)
    return e1, np.cross(nu, e1)


def calculadora_recta_plano(modo, A_txt, B_txt, C_txt, n_txt, d_txt, P0_txt, u_txt, Q_txt, n_cifras):
    try:
        if modo == "plano por tres puntos":
            pts = [vector(t, 3, e) for t, e in ((A_txt, "A"), (B_txt, "B"), (C_txt, "C"))]
            n, d = plano_tres_puntos(*pts)
        else:
            pts = []
            n, d = _normal(n_txt), numero(d_txt, "d")
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"Plano: {plano_texto(n, d, n_cifras)}      (normal n = {fmt(n, n_cifras)}, d = {_num_txt(d, n_cifras)})")
    if pts:
        print("   n = AB × AC simplificado. Comprobación: n · A, n · B, n · C =",
              ", ".join(_num_txt(_punto(n, P), n_cifras) for P in pts), " (las tres deben dar d)")
        print("   (en este modo no se usan los campos n y d)")
    else:
        print("   (en este modo no se usan los campos A, B y C)")

    Q = H = P0 = X = None
    try:
        Q = vector(Q_txt, 3, "Q")
        D, lado = distancia(Q, n, d)
        print(f"Distancia de Q = {fmt(Q, n_cifras)} al plano: D = {_num_txt(D, n_cifras)}; Q queda {lado}.")
        H = _pie(Q, n, d)                               # pie de la perpendicular, para el dibujo
    except ValueError as err:
        print("Revisa la entrada:", err); Q = None
    try:
        P0, u = vector(P0_txt, 3, "P0"), vector(u_txt, 3, "la dirección u")
        tipo, t, X = interseccion(P0, u, n, d)
    except ValueError as err:
        print("Revisa la entrada:", err); P0 = None
    if P0 is not None:
        print(f"Recta: r(t) = P0 + t u = {fmt(P0, n_cifras)} + t {fmt(u, n_cifras)}")
        print("   paramétricas:", _parametricas(P0, u, n_cifras))
        if tipo == "punto":
            print(f"La recta corta al plano en t* = {_num_txt(t, n_cifras)}, en el punto {_fmt_txt(X, n_cifras)}.")
            if t < 0:
                print("t* < 0: el corte queda antes de P0 (detrás, si la recta es un rayo)")
        elif tipo == "paralela":
            print("La recta es paralela al plano y no lo toca: n · u = 0 y n · P0 ≠ d.")
        else:
            print("La recta está contenida en el plano: n · u = 0 y n · P0 = d.")

    # Figura: un trozo del plano centrado en lo que interesa, la recta, el corte y la distancia
    interes = pts + [p for p in (Q, H, P0, X) if p is not None]
    if not interes:
        interes = [_pie(np.zeros(3), n, d)]
    centro = _pie(np.mean(interes, axis=0), n, d)
    R = 1.2 * max(_norma(p - centro) for p in interes) or 1.0
    e1, e2 = _base_plano(n)
    s1, s2 = np.meshgrid([-R, R], [-R, R])
    malla = [centro[i] + s1 * e1[i] + s2 * e2[i] for i in range(3)]
    fig = plt.figure(figsize=(6.5, 6))
    ax = fig.add_subplot(projection="3d")
    ax.plot_surface(*malla, color="tab:blue", alpha=0.25, shade=False)
    ax.quiver(*centro, *(0.5 * R * n / _norma(n)), color="tab:blue", arrow_length_ratio=0.15, label="dirección de n")
    dibujados = [centro, *(np.array(malla).reshape(3, -1).T)]
    for P, e in zip(pts, "ABC"):
        ax.scatter(*P, color="tab:blue", s=15); ax.text(*P, f"  {e}")
    if P0 is not None:
        L = R / _norma(u)
        t0, t1 = (min(0.0, t), max(0.0, t)) if tipo == "punto" else (0.0, 0.0)
        extremos = np.array([P0 + (t0 - L) * u, P0 + (t1 + L) * u])
        ax.plot(*extremos.T, color="black", lw=1.5, label="recta P0 + t u")
        mismo = Q is not None and np.array_equal(P0, Q)
        ax.scatter(*P0, color="black", s=20); ax.text(*P0, "  P0 = Q" if mismo else "  P0")
        dibujados += list(extremos)
        if tipo == "punto":
            ax.scatter(*X, color="red", s=45, label=f"corte (t* = {_num_txt(t, 3)})")
    if Q is not None:
        ax.scatter(*Q, color="green", s=20)
        if P0 is None or not np.array_equal(P0, Q):
            ax.text(*Q, "  Q")
        ax.plot(*np.array([Q, H]).T, ":", color="green", lw=1.5, label=f"distancia D = {_num_txt(D, 3)}")
    _ejes_iguales(ax, dibujados + interes)
    _vista(ax, [u] if P0 is not None else [], normal=n)   # el plano no se ve de canto ni de frente
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.set_title(f"Plano {plano_texto(n, d, 3)}")
    ax.legend(loc="upper left", fontsize=8)
    plt.show()


_ancho = {"description_width": "initial"}
widgets.interact(calculadora_recta_plano,
    modo=widgets.Dropdown(options=["plano por tres puntos", "plano por normal y d"], value="plano por tres puntos",
                          description="modo"),
    A_txt=widgets.Text(value="1, 2, 0", description="A =", continuous_update=False),
    B_txt=widgets.Text(value="3, 0, 1", description="B =", continuous_update=False),
    C_txt=widgets.Text(value="0, 1, 2", description="C =", continuous_update=False),
    n_txt=widgets.Text(value="2, 1, 2", description="n =", continuous_update=False),
    d_txt=widgets.Text(value="12", description="d =", continuous_update=False),
    P0_txt=widgets.Text(value="0, 0, 0", description="P0 (recta) =", continuous_update=False, style=_ancho),
    u_txt=widgets.Text(value="1, 1, 1", description="u (recta) =", continuous_update=False, style=_ancho),
    Q_txt=widgets.Text(value="2, 3, 4", description="Q (distancia) =", continuous_update=False, style=_ancho),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras sig.", continuous_update=False));

# %%
# Casos de prueba de la sección 5 (resultado conocido)
def _plano_libro():
    """El plano por (1, 2, 0), (3, 0, 1) y (0, 1, 2) tiene n ∥ (3, 5, 4) y pasa por los tres puntos."""
    puntos = ((1, 2, 0), (3, 0, 1), (0, 1, 2))
    n, d = plano_tres_puntos(*puntos)
    return (paralelos(n, (3, 5, 4)) and all(cerca(punto(n, P), d) for P in puntos)
            and plano_texto(n, d) == "3x + 5y + 4z = 13"
            and misma_ecuacion((6, 10, 8, 26), (3, 5, 4, 13)) and not misma_ecuacion((3, 5, 4, 12), (3, 5, 4, 13)))


PRUEBAS_5 = [
    ("plano por (1, 2, 0), (3, 0, 1), (0, 1, 2): 3x + 5y + 4z = 13, pasa por los tres; (6, 10, 8, 26) es el mismo",
     _plano_libro),
    ("distancia de (4, 5, 6) a x + 2y + 2z = 12 = 14/3, del lado hacia donde apunta n",
     lambda: cerca(distancia((4, 5, 6), (1, 2, 2), 12)[0], 14 / 3)
             and distancia((4, 5, 6), (1, 2, 2), 12)[1] == "del lado hacia donde apunta n"),
    ("recta (0, 0, 1) + t (1, 1, 1) y plano 2x + y + 2z = 12: corta en t = 2, en (2, 2, 3)",
     lambda: interseccion((0, 0, 1), (1, 1, 1), (2, 1, 2), 12)[0] == "punto"
             and cerca(interseccion((0, 0, 1), (1, 1, 1), (2, 1, 2), 12)[1], 2)
             and np.allclose(interseccion((0, 0, 1), (1, 1, 1), (2, 1, 2), 12)[2], (2, 2, 3))),
    ("con u = (1, 0, -1) la recta es paralela al plano",
     lambda: interseccion((0, 0, 1), (1, 0, -1), (2, 1, 2), 12)[0] == "paralela"),
    ("los puntos alineados (1, 2, 3), (3, 5, 4), (5, 8, 5) se rechazan",
     lambda: _rechaza(lambda: plano_tres_puntos((1, 2, 3), (3, 5, 4), (5, 8, 5)))),
    ("recta (1, 1, 0) + t (1, 0, 0) contenida en el plano z = 0",
     lambda: interseccion((1, 1, 0), (1, 0, 0), (0, 0, 1), 0)[0] == "contenida"),
]
for nombre, prueba in PRUEBAS_5:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 5).** Encuentra a mano la ecuación $ax+by+cz=d$ del plano que pasa por $(2,0,0)$, $(0,4,0)$ y $(0,0,4)$. Escribe los cuatro números $a, b, c, d$. Tu ecuación puede diferir de la de la calculadora en un factor: cualquier múltiplo no nulo describe el mismo plano.

# %%
mi_plano = None     # escribe a, b, c, d de tu ecuación ax + by + cz = d, por ejemplo: "1, -1, 2, 5"

if mi_plano is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    n_ref, d_ref = plano_tres_puntos((2, 0, 0), (0, 4, 0), (0, 0, 4))
    print("La calculadora da:", plano_texto(n_ref, d_ref))
    try:
        ok = misma_ecuacion(mi_plano, (*n_ref, d_ref), tol=1e-3)
    except ValueError as err:
        if isinstance(mi_plano, str) and re.search(r"[xyz=]", mi_plano):
            err = "escribe solo los cuatro números a, b, c, d; por ejemplo, para 2x + 3y - z = 5 escribe \"2, 3, -1, 5\"."
        print("Revisa tu valor:", err); ok = None
    if ok is not None:
        print("coinciden: describen el mismo plano" if ok else
              "NO coinciden; usa n = AB × AC, d = n · A, y comprueba que los tres puntos cumplan tu ecuación")
