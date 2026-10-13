# ID: ALG-U4-NB03
# Notebook: alg/u4_espacios.ipynb · sección 4.3 base, dimensión, coordenadas y cambio de base
# Repositorio: alg/u4_espacios/03_base_coordenadas.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Base, dimensión, coordenadas y cambio de base
#
# **Base.** Una *base* de un subespacio $H$ es un conjunto de vectores de $H$ que es **linealmente independiente** y **genera** $H$. Todas las bases de $H$ tienen el mismo número de vectores: la *dimensión* de $H$, $\dim H$ ($\dim\{\mathbf 0\}=0$).
#
# - La **base estándar** de $\mathbb R^n$ es $\mathbf e_1=(1,0,\dots,0),\ \dots,\ \mathbf e_n=(0,\dots,0,1)$, así que $\dim\mathbb R^n=n$.
# - En $\mathbb R^n$, **$n$ vectores independientes ya son base** (generan todo $\mathbb R^n$); con $n$ vectores, eso equivale a $\det[\,\mathbf b_1\ \cdots\ \mathbf b_n\,]\neq 0$.
#
# **Los tres subespacios de una matriz $A$ de $m\times n$.** Al escalonar $A$:
#
# - **base de $\operatorname{Col}(A)$**: las columnas pivote **de la $A$ original** (no las de la forma reducida: escalonar cambia el espacio columna);
# - **base de $\operatorname{Nul}(A)$**: los vectores de la solución general de $A\mathbf x=\mathbf 0$, uno por variable libre;
# - **base del espacio renglón $\operatorname{Ren}(A)$**: los renglones no nulos de la forma escalonada (las operaciones de renglón no cambian el espacio renglón).
#
# $$\operatorname{rango}(A)=\dim\operatorname{Col}(A)=\dim\operatorname{Ren}(A),\qquad \operatorname{rango}(A)+\dim\operatorname{Nul}(A)=n\quad\text{(teorema de la dimensión)}.$$
#
# **Coordenadas.** Si $B=\{\mathbf b_1,\dots,\mathbf b_n\}$ es base, todo $\mathbf x$ se escribe **de una sola manera** como $\mathbf x=c_1\mathbf b_1+\dots+c_n\mathbf b_n$ (de dos maneras, la resta daría una relación de dependencia). Los coeficientes son las *coordenadas de $\mathbf x$ en la base $B$*, $[\mathbf x]_B=(c_1,\dots,c_n)$. Con $P_B=[\,\mathbf b_1\ \cdots\ \mathbf b_n\,]$:
#
# $$\mathbf x=P_B\,[\mathbf x]_B,\qquad [\mathbf x]_B=P_B^{-1}\,\mathbf x.$$
#
# **Cambio de base.** Si $B$ y $C$ son bases de $\mathbb R^n$, de $P_B[\mathbf x]_B=\mathbf x=P_C[\mathbf x]_C$ sale
#
# $$[\mathbf x]_C=P_{C\leftarrow B}\,[\mathbf x]_B,\qquad P_{C\leftarrow B}=P_C^{-1}P_B .$$
#
# Las columnas de $P_{C\leftarrow B}$ son las coordenadas en $C$ de los vectores de $B$.
#
# **Base de un subespacio.** Si $B$ tiene $k<n$ vectores independientes, es base del subespacio que generan (un plano en $\mathbb R^3$, por ejemplo). Un $\mathbf x$ de ese subespacio tiene $k$ coordenadas; un $\mathbf x$ fuera de él no tiene coordenadas en $B$.

# %%
def _revisa_base(P, nombre="la base"):
    """ValueError si las columnas de P (los vectores de la base) no son independientes."""
    n, k = P.shape
    if k > n:
        raise ValueError(f"no es base: {nombre} tiene {k} vectores en {_R(n)}, y más de {n} vectores en {_R(n)} "
                         "siempre son dependientes")
    r = _elimina(P, k, con_b=False)["r"]
    if r < k:
        detalle = f" (det = {P.det()})" if k == n else f" (rango {r} < {k} vectores)"
        raise ValueError(f"no es base: los vectores son dependientes{detalle}")


def coordenadas(base, x):
    """[x]_B exacto: los coeficientes c con x = c1·b1 + ... + ck·bk (B escrita como lista de vectores, uno por
    renglón). Si B tiene n vectores dependientes en Rⁿ: ValueError 'no es base: los vectores son dependientes'.
    Si B tiene k < n vectores independientes (base de un subespacio), da las coordenadas si x está en el
    subespacio; si no, ValueError 'x no está en el subespacio generado por la base'."""
    P = lista_vectores(base, "la base")
    x = vector(x, "x")
    n, k = P.shape
    if x.rows != n:
        raise ValueError(f"x tiene {_componentes(x.rows)} y los vectores de la base tienen {n}: deben coincidir")
    _revisa_base(P)
    el = _elimina(P.row_join(x), k, con_b=True)
    if el["imposible"] is not None:
        raise ValueError("x no está en el subespacio generado por la base")
    return el["R"][:k, k]


def desde_coordenadas(base, c):
    """x = P_B·[x]_B = c1·b1 + ... + ck·bk."""
    P = lista_vectores(base, "la base")
    c = vector(c, "las coordenadas")
    if c.rows != P.cols:
        raise ValueError(f"escribiste {_componentes(c.rows).replace('componente', 'coordenada')} y la base tiene "
                         f"{P.cols} {'vector' if P.cols == 1 else 'vectores'}: hace falta una coordenada por vector")
    _revisa_base(P)
    return P * c


def matriz_cambio(B, C):
    """P_{C←B} = P_C⁻¹·P_B, exacta (B y C, bases de Rⁿ escritas como listas de vectores)."""
    PB, PC = lista_vectores(B, "la base B"), lista_vectores(C, "la base C")
    if PB.rows != PC.rows:
        raise ValueError(f"los vectores de B tienen {PB.rows} componentes y los de C {PC.rows}: deben ser del mismo espacio")
    for P, letra in ((PB, "B"), (PC, "C")):
        n, k = P.shape
        if k != n:
            raise ValueError(f"la base {letra} tiene {k} {'vector' if k == 1 else 'vectores'} en {_R(n)}: para cambiar de "
                             f"base en {_R(n)}, cada base necesita {n} vectores")
        try:
            _revisa_base(P, f"la base {letra}")
        except ValueError as err:
            raise ValueError(f"{letra} {err}") from None
    return PC.inv() * PB


def base_col_nul(A):
    """Bases de Col(A), Nul(A) y Ren(A) de una matriz A escrita POR RENGLONES (hasta 6×6). Diccionario con reducida,
    escalonada, pivote (columnas pivote, desde 1), base_col (columnas pivote de la A ORIGINAL), base_nul (vectores
    exactos de la solución general de A x = 0), base_ren (renglones no nulos de la forma escalonada, como vectores),
    rango, nulidad, n y comprobacion (rango + nulidad == n)."""
    A = matriz_limitada(A, "A")
    m, n = A.shape
    el = _elimina(A, n, con_b=False)
    piv = el["pivotes"]
    base_nul = _direcciones(el["R"], piv, n)
    return {"A": A, "escalonada": el["E"], "reducida": el["R"], "pivote": [c + 1 for c in piv],
            "base_col": [A[:, c] for c in piv], "base_nul": base_nul,
            "base_ren": [el["E"].row(i).T for i in range(el["r"])], "rango": el["r"], "nulidad": len(base_nul),
            "n": n, "comprobacion": el["r"] + len(base_nul) == n}


def _informe_coordenadas(P, x, c, n_cifras, modo="x → coordenadas", letra="b"):
    """Informe de las coordenadas; letra es el nombre de los vectores de la base (b1, b2... o a1, a2...)."""
    n, k = P.shape
    nombres = [f"{letra}{j}" for j in range(1, k + 1)]
    print("Base B: " + ", ".join(f"{nm} = {_tupla(P[:, j])}" for j, nm in enumerate(nombres)))
    muestra("P_B", P, n_cifras)
    if k == n:
        d = P.det()
        igual = _igual("det P_B", d, n_cifras)
        print(igual + (" ≠ 0: " if "(" not in igual else ", distinto de 0: ") +
              (f"el vector es independiente, así que B es base de {_R(n)}." if n == 1 else
              f"los {n} vectores son independientes, así que B es base de {_R(n)}."))
    else:
        print(f"B tiene {k} {'vector' if k == 1 else 'vectores'} {'independiente' if k == 1 else 'independientes'} en "
              f"{_R(n)}: es base de {_subespacio(k, n)}, no de todo {_R(n)}.")
    if modo == "x → coordenadas":
        if k == n:
            muestra("P_B⁻¹", P.inv(), n_cifras)
            print(f"[x]_B = P_B⁻¹·x = {_tupla(c, n_cifras)}")
        else:
            print(f"[x]_B sale de resolver P_B c = x (los rangos coinciden: x está en el subespacio): {_tupla(c, n_cifras)}")
    else:
        print(f"x = P_B·[x]_B = {_tupla(x, n_cifras)}")
    print(f"Comprobación: {_texto_combinacion(list(c), nombres=nombres)} = {_tupla(P * c)}" +
          (" = x" if P * c == x else "  ≠ x: ¡revisa!"))


_TOPE_DIBUJO = 10 ** 6          # más allá de esto las rejillas ya no se pueden leer: no se dibujan


def _paso_rejilla(lo, hi, maximo=40):
    """Paso entero para trazar a lo más unas 40 líneas entre lo y hi."""
    return max(1, math.ceil((hi - lo) / maximo))


def dibuja_rejillas(base, x, letra="b"):
    """Solo R²: rejilla estándar (gris claro, continua) y rejilla de la base B (discontinua), los vectores de la base
    (el primero continuo, el segundo discontinuo), el punto x y el camino c1·b1 seguido de c2·b2 (punteado grueso);
    el título da las coordenadas en las dos bases. Cada rejilla lleva a lo más unas 40 líneas; si las coordenadas
    pasan de 10⁶, no dibuja y lo dice."""
    P = lista_vectores(base, "la base")
    xv = vector(x, "x")
    if P.shape != (2, 2):
        raise ValueError(f"las rejillas se dibujan para una base de R² con 2 vectores (esta tiene {P.cols} "
                         f"{'vector' if P.cols == 1 else 'vectores'} de {P.rows} componentes)")
    c = coordenadas(P, xv)
    b1, b2 = _flotante(P[:, 0]).ravel(), _flotante(P[:, 1]).ravel()
    xf, cf = _flotante(xv).ravel(), _flotante(c).ravel()
    puntos = np.vstack([np.zeros(2), b1, b2, xf, cf[0] * b1])
    L = 1.25 * max(1.0, float(np.abs(puntos).max())) + 0.5
    if not math.isfinite(L) or L > _TOPE_DIBUJO or float(np.abs(cf).max()) > _TOPE_DIBUJO:
        print(f"No se dibujan las rejillas: las coordenadas pasan de {_TOPE_DIBUJO:,} y la figura no se podría leer.")
        return c
    n1, n2 = f"{letra}1", f"{letra}2"
    fig, ax = plt.subplots(figsize=(6.4, 6.4))
    K = math.ceil(L)
    paso = _paso_rejilla(-K, K)
    for k in range(-K - (-K % paso), K + 1, paso):          # rejilla estándar (a lo más unas 40 líneas)
        ax.axvline(k, color="0.82", lw=0.8, zorder=0, label=None)
        ax.axhline(k, color="0.82", lw=0.8, zorder=0)
    ax.plot([], [], "-", color="0.82", lw=0.8, label="rejilla estándar" + (f" (cada {paso})" if paso > 1 else ""))
    esquinas = np.array([[-L, -L], [L, -L], [L, L], [-L, L]])
    en_B = (_flotante(P.inv()) @ esquinas.T).T                  # esquinas de la ventana en coordenadas de B (P_B⁻¹ exacta)
    pasos_B = []
    for eje, (d, otro) in enumerate(((b2, b1), (b1, b2))):      # líneas i·b1 + t·b2 y j·b2 + t·b1
        lo, hi = math.floor(en_B[:, eje].min()), math.ceil(en_B[:, eje].max())
        paso_B = _paso_rejilla(lo, hi)
        pasos_B.append(paso_B)
        t = np.array([-1.0, 1.0]) * 4 * L / max(np.hypot(*d), 1e-12)
        for i in range(lo - (lo % paso_B), hi + 1, paso_B):
            base_pt = i * otro
            ax.plot(base_pt[0] + t * d[0], base_pt[1] + t * d[1], "--", color=_COLORES[0], lw=0.9, alpha=0.8)
    ax.plot([], [], "--", color=_COLORES[0], lw=0.9, label="rejilla de B" + (" (no todas las líneas)" if max(pasos_B) > 1 else ""))
    _flecha(ax, (0, 0), b1, _COLORES[1], "-", 2.6, n1, (6, -14))
    _flecha(ax, (0, 0), b2, _COLORES[2], "--", 2.6, n2)
    ax.plot([], [], "-", color=_COLORES[1], lw=2.6, label=f"{n1} (flecha continua)")
    ax.plot([], [], "--", color=_COLORES[2], lw=2.6, label=f"{n2} (flecha discontinua)")
    p1 = cf[0] * b1
    camino = f"camino: {_termino(c[0], n1, True)} y luego {_termino(c[1], n2, True)}"
    ax.plot(*np.array([[0, 0], p1]).T, ":", color="black", lw=3, label=camino)
    ax.plot(*np.array([p1, xf]).T, ":", color="black", lw=3)
    ax.plot(*p1, "s", color="black", ms=6)
    ax.plot(*xf, "o", color="black", ms=9, zorder=5, label=f"x = {_tupla(xv)}")
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect("equal")
    ax.axhline(0, color="black", lw=0.8); ax.axvline(0, color="black", lw=0.8)
    ax.set_xlabel("x (estándar)"); ax.set_ylabel("y (estándar)")
    ax.set_title(textwrap.fill(f"En la base estándar: x = {_tupla(xv)};  en la base B: [x]_B = {_tupla(c)}", 70), fontsize=10)
    ax.legend(loc="lower right", fontsize=8)
    plt.show()
    return c

# %% [markdown]
# **Recorrido paso a paso.** La semilla del libro: en una ciudad las calles van **al este**, $\mathbf a_1=(1,0)$, y **al noreste**, $\mathbf a_2=(1,1)$; la base es $B=\{\mathbf a_1,\mathbf a_2\}$. Para llegar a $\mathbf x=(3,2)$ (3 al este y 2 al norte en el plano usual) caminas $c_1$ cuadras al este y $c_2$ al noreste: $c_1(1,0)+c_2(1,1)=(3,2)$, así que $c_2=2$ y $c_1=1$, es decir, $[\mathbf x]_B=(1,2)$. En la gráfica, la rejilla de $B$ (discontinua) es el plano de calles: $\mathbf x$ queda en la esquina 1 al este y 2 al noreste.

# %%
_PB = lista_vectores("1, 0; 1, 1", "la base")
_xB = coordenadas(_PB, "3, 2")
_informe_coordenadas(_PB, vector("3, 2"), _xB, 4, letra="a")
dibuja_rejillas("1, 0; 1, 1", "3, 2", letra="a");

# %% [markdown]
# **Calculadora de coordenadas.** Escribe los vectores de la base $B$ **uno por renglón** y un vector. En el modo «x → coordenadas», el vector es $\mathbf x$ y la calculadora da $[\mathbf x]_B$; en el modo «coordenadas → x», el vector es $[\mathbf x]_B$ y la calculadora da $\mathbf x$. Con una base de $\mathbb R^2$ dibuja las dos rejillas. Prueba una base que no lo sea (`1, 2; 2, 4`) y una base de un plano de $\mathbb R^3$ con un $\mathbf x$ fuera de él.

# %%
def calculadora_coordenadas(base_txt, vector_txt, modo, dibujar, n_cifras):
    with bloque():                 # el texto aparece completo, de una sola vez; el dibujo va después
        try:
            P = lista_vectores(base_txt, "la base")
            if modo == "x → coordenadas":
                x = vector(vector_txt, "x")
                c = coordenadas(P, x)
            else:
                c = vector(vector_txt, "[x]_B")
                x = desde_coordenadas(P, c)
        except ValueError as err:
            print("Revisa la entrada:", err); return
        _informe_coordenadas(P, x, c, n_cifras, modo)
    if dibujar:
        if P.shape == (2, 2):
            dibuja_rejillas(P, x)
        else:
            print("Las rejillas se dibujan solo para una base de R² (dos vectores de 2 componentes).")


_CONTROLES_3A = dict(
    base_txt=widgets.Textarea(value="1, 0\n1, 1", description="base B:", continuous_update=False,
                              layout=widgets.Layout(width="300px", height="100px")),
    vector_txt=widgets.Text(value="3, 2", description="x o [x]_B:", continuous_update=False,
                            style={"description_width": "initial"}),
    modo=widgets.Dropdown(options=["x → coordenadas", "coordenadas → x"], value="x → coordenadas", description="modo:"),
    dibujar=widgets.Checkbox(value=True, description="dibujar rejillas (R²)"),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_coordenadas, **_CONTROLES_3A);

# %% [markdown]
# **Calculadora de los subespacios de una matriz.** Escribe $A$ **por renglones** (como en la unidad 3). La calculadora da las bases de $\operatorname{Col}(A)$, $\operatorname{Nul}(A)$ y $\operatorname{Ren}(A)$, el rango, la nulidad y comprueba el teorema de la dimensión.

# %%
def calculadora_col_nul(A_txt, n_cifras):
    with bloque():
        try:
            res = base_col_nul(A_txt)
        except ValueError as err:
            print("Revisa la entrada:", err); return
        A, n = res["A"], res["n"]
        muestra("A", A, n_cifras)
        muestra("forma escalonada", res["escalonada"], n_cifras)
        muestra("forma escalonada reducida", res["reducida"], n_cifras)
        piv = res["pivote"]
        print(f"{'Columna pivote' if len(piv) == 1 else 'Columnas pivote'}: {_lista(piv) if piv else 'ninguna'}")
        if piv:
            print(f"Base de Col(A) ({'columna' if len(piv) == 1 else 'columnas'} {_lista(piv)} de la A ORIGINAL): " +
                  ", ".join(_tupla(v, n_cifras) for v in res["base_col"]))
        else:
            print("Base de Col(A): ninguna, porque A no tiene columnas pivote: Col(A) = {0}")
        print("Base de Nul(A) (solución general de A x = 0, un vector por variable libre): " +
              (", ".join(_tupla(v, n_cifras) for v in res["base_nul"]) if res["base_nul"] else "ninguna: Nul(A) = {0}"))
        print("Base de Ren(A) (renglones no nulos de la forma escalonada): " +
              (", ".join(_tupla(v, n_cifras) for v in res["base_ren"]) if res["base_ren"] else "ninguna: Ren(A) = {0}"))
        print(f"rango(A) = dim Col(A) = dim Ren(A) = {res['rango']};  nulidad = dim Nul(A) = {res['nulidad']}")
        print(f"Teorema de la dimensión: rango + nulidad = {res['rango']} + {res['nulidad']} = "
              f"{res['rango'] + res['nulidad']} = n (columnas de A) = {n}" + ("  ✓" if res["comprobacion"] else "  ✗"))
        ceros = all(A * v == sp.zeros(A.rows, 1) for v in res["base_nul"])
        if res["base_nul"]:
            print("Comprobación: A·v = 0 para cada vector de la base de Nul(A)." if ceros else "Comprobación FALLÓ: A·v ≠ 0.")


_CONTROLES_3B = dict(
    A_txt=widgets.Textarea(value="1, 2, 3; 2, 4, 7; 3, 6, 10", description="A =", continuous_update=False,
                           layout=widgets.Layout(width="340px", height="110px")),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_col_nul, **_CONTROLES_3B);

# %% [markdown]
# **Calculadora de cambio de base.** Escribe las bases $B$ y $C$ de $\mathbb R^n$ (los vectores separados por `;`) y unas coordenadas $[\mathbf x]_B$. La calculadora arma $P_{C\leftarrow B}=P_C^{-1}P_B$ y da $[\mathbf x]_C$.

# %%
def calculadora_cambio(B_txt, C_txt, xB_txt, n_cifras):
    with bloque():
        try:
            M = matriz_cambio(B_txt, C_txt)
            PB, PC = lista_vectores(B_txt, "la base B"), lista_vectores(C_txt, "la base C")
            xB = vector(xB_txt, "[x]_B")
            if xB.rows != PB.cols:
                raise ValueError(f"[x]_B tiene {_componentes(xB.rows)} y la base B tiene {PB.cols} vectores: "
                                 "hace falta una coordenada por vector")
        except ValueError as err:
            print("Revisa la entrada:", err); return
        muestra("P_B", PB, n_cifras)
        muestra("P_C", PC, n_cifras)
        muestra("P_C⁻¹", PC.inv(), n_cifras)
        muestra("P_{C←B} = P_C⁻¹·P_B", M, n_cifras)
        xC = M * xB
        x = PB * xB
        print(f"[x]_C = P_{{C←B}}·[x]_B = {_tupla(xC, n_cifras)}")
        print(f"Comprobación: x = P_B·[x]_B = {_tupla(x)} y P_C·[x]_C = {_tupla(PC * xC)}" +
              ("  (el mismo vector)" if PC * xC == x else "  ¡no coinciden!"))


_CONTROLES_3C = dict(
    B_txt=widgets.Text(value="1, 0; 1, 1", description="B:", continuous_update=False),
    C_txt=widgets.Text(value="1, 2; -1, 1", description="C:", continuous_update=False),
    xB_txt=widgets.Text(value="1, 2", description="[x]_B =", continuous_update=False),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_cambio, **_CONTROLES_3C);

# %%
# Casos de prueba de la sección 3 (resultado conocido)
PRUEBAS_3 = [
    ("base (1, 0), (1, 1), x = (3, 2) → [x]_B = (1, 2)",
     lambda: coordenadas("1,0; 1,1", "3,2") == vector("1, 2")),
    ("base (1, 2), (3, 4), x = (5, 6) → [x]_B = (-1, 2)",
     lambda: coordenadas("1,2; 3,4", "5,6") == vector("-1, 2")),
    ("(1, 2), (2, 4) no es base → ValueError «no es base: los vectores son dependientes»",
     lambda: _mensaje(lambda: coordenadas("1,2; 2,4", "1,1")).startswith("no es base: los vectores son dependientes")),
    ("P_(C←B) con B = (1, 0), (1, 1) y C = (1, 1), (1, -1) → [[1/2, 1], [1/2, 0]]; aplicada a (1, 2) da (5/2, 1/2)",
     lambda: (lambda M: M == sp.Matrix([[sp.Rational(1, 2), 1], [sp.Rational(1, 2), 0]])
                        and M * vector("1, 2") == vector("5/2, 1/2"))(matriz_cambio("1,0; 1,1", "1,1; 1,-1"))),
    ("base_col_nul(1, 2, 3; 2, 4, 7; 3, 6, 10) → pivote [1, 3], Col: (1, 2, 3), (3, 7, 10), Nul: (-2, 1, 0), 2 + 1 = 3",
     lambda: (lambda r: r["pivote"] == [1, 3] and r["base_col"] == [vector("1, 2, 3"), vector("3, 7, 10")]
                        and len(r["base_nul"]) == 1 and sp.Matrix.hstack(r["base_nul"][0], vector("-2, 1, 0")).rank() == 1
                        and (r["rango"], r["nulidad"], r["n"]) == (2, 1, 3) and r["comprobacion"])(
         base_col_nul("1,2,3; 2,4,7; 3,6,10"))),
    ("base de un plano (2, 1, 0), (-1, 0, 1): x = (3, 2, 1) → (2, 1); x = (1, 0, 0) → no está en el subespacio",
     lambda: coordenadas("2,1,0; -1,0,1", "3,2,1") == vector("2, 1")
             and _mensaje(lambda: coordenadas("2,1,0; -1,0,1", "1,0,0")) == "x no está en el subespacio generado por la base"),
    ("base (1/2, 1/2), (1/2, -1/2), x = (30, 10) → [x]_B = (40, 20)",
     lambda: coordenadas("1/2,1/2; 1/2,-1/2", "30,10") == vector("40, 20")),
    ("ida y vuelta con la base (1, 2), (3, 4) y x = (7, -3): desde_coordenadas(B, coordenadas(B, x)) = x",
     lambda: desde_coordenadas("1,2; 3,4", coordenadas("1,2; 3,4", "7,-3")) == vector("7, -3")),
    ("base_col_nul de la identidad 3×3 → rango 3, base de Nul vacía",
     lambda: (lambda r: r["rango"] == 3 and r["base_nul"] == [] and r["comprobacion"])(base_col_nul(sp.eye(3)))),
    ("base de Ren(1, 2, 3; 2, 4, 7; 3, 6, 10) → los dos renglones no nulos de la escalonada: (1, 2, 3) y (0, 0, 1)",
     lambda: base_col_nul("1,2,3; 2,4,7; 3,6,10")["base_ren"] == [vector("1, 2, 3"), vector("0, 0, 1")]),
    ("desde_coordenadas con una base dependiente, (1, 2), (2, 4) → ValueError «no es base»",
     lambda: _mensaje(lambda: desde_coordenadas("1,2; 2,4", "1, 1")).startswith("no es base")),
    ("prueba de humo: la calculadora de coordenadas corre con sus valores iniciales, en los dos modos",
     lambda: _humo(calculadora_coordenadas, _CONTROLES_3A)
             and _humo(calculadora_coordenadas, _CONTROLES_3A, vector_txt="1, 2", modo="coordenadas → x")),
    ("prueba de humo: la calculadora de Col, Nul y Ren corre con sus valores iniciales",
     lambda: _humo(calculadora_col_nul, _CONTROLES_3B)),
    ("prueba de humo: la calculadora de cambio de base corre con sus valores iniciales",
     lambda: _humo(calculadora_cambio, _CONTROLES_3C)),
    ("rejillas con coordenadas enormes: no se cuelga (x = (10000, 1) con b2 = (1, 0.0001); x = (1e15, 0))",
     lambda: (lambda t0: _humo(calculadora_coordenadas, _CONTROLES_3A, base_txt="1, 0; 1, 0.0001", vector_txt="10000, 1")
                         and "No se dibujan" in _captura(calculadora_coordenadas, "1, 0; 1, 1", "1e15, 0",
                                                         "x → coordenadas", True, 4)[0]
                         and time.time() - t0 < 20)(time.time())),
]
corre_pruebas(PRUEBAS_3, 3)

# %% [markdown]
# **Contrasta (sección 3).** Calcula a mano las coordenadas de $\mathbf x=(4,1)$ en la base $B=\{(1,1),(1,-1)\}$: plantea $c_1(1,1)+c_2(1,-1)=(4,1)$ y resuelve. Escribe `mis_coordenadas = (c1, c2)`; puedes usar fracciones como texto (`"7/2"`) o decimales.

# %%
mis_coordenadas = None     # dos números, por ejemplo (3, -1/2) escrito como ("3", "-1/2") o (3, -0.5)

if mis_coordenadas is None:     # la referencia se muestra después de tu cálculo
    print("Falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda.")
else:
    try:
        mias = vector(mis_coordenadas, "tus coordenadas")
    except ValueError as err:
        print(f"Revisa tus coordenadas: {str(err).rstrip('.')}.")
    else:
        if mias.rows != 2:          # respuesta imposible: la referencia todavía no se muestra
            print(f"La base tiene dos vectores, así que son dos coordenadas (escribiste {mias.rows}): revisa tu "
                  "respuesta y vuelve a ejecutar la celda.")
        else:
            ref = coordenadas("1, 1; 1, -1", "4, 1")
            print(f"La calculadora da: [x]_B = {_tupla(ref, 4)}.")
            if mias == ref:
                print("Coinciden.")
            elif all(abs(float(a - b)) <= 0.005 for a, b in zip(mias, ref)):
                print("Coinciden salvo redondeo (escribe la fracción exacta si puedes).")
            else:
                print("NO coinciden: suma y resta las dos ecuaciones, c1 + c2 = 4 y c1 - c2 = 1; comprueba que "
                      "c1·(1, 1) + c2·(1, -1) dé (4, 1).")
            print("\nEl procedimiento completo:")
            _P = lista_vectores("1, 1; 1, -1", "la base")
            _informe_coordenadas(_P, vector("4, 1"), ref, 4)
