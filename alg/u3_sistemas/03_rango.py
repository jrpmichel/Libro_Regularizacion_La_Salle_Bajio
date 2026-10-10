# ID: ALG-U3-NB03
# Notebook: alg/u3_sistemas.ipynb · sección 3.3 ecuaciones vectorial y matricial; rango
# Repositorio: alg/u3_sistemas/03_rango.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Ecuaciones vectorial y matricial; rango
#
# El mismo sistema se escribe de tres maneras. Con $\mathbf a_j$ la columna $j$ de $A$:
#
# $$\text{ecuaciones: }\ a_{i1}x_1+\dots+a_{in}x_n=b_i\ \ (i=1,\dots,m),\qquad \text{matricial: }\ A\mathbf x=\mathbf b,\qquad \text{vectorial: }\ x_1\mathbf a_1+x_2\mathbf a_2+\dots+x_n\mathbf a_n=\mathbf b.$$
#
# - **Lectura por renglones.** Cada ecuación es una recta (o un plano) y la solución es el punto común a todas.
# - **Lectura por columnas.** El sistema pregunta con qué coeficientes hay que combinar las columnas de $A$ para llegar a $\mathbf b$. Una suma $c_1\mathbf a_1+\dots+c_n\mathbf a_n$ es una **combinación lineal** de las columnas.
#
# Por ejemplo, $\begin{bmatrix}2&1\\1&3\end{bmatrix}\begin{bmatrix}2\\1\end{bmatrix}=2\begin{bmatrix}2\\1\end{bmatrix}+1\begin{bmatrix}1\\3\end{bmatrix}=\begin{bmatrix}5\\5\end{bmatrix}$: el producto $A\mathbf x$ **es** la combinación de las columnas con los coeficientes de $\mathbf x$. Por eso:
#
# > $A\mathbf x=\mathbf b$ tiene solución si y solo si $\mathbf b$ es combinación lineal de las columnas de $A$, y los coeficientes de la combinación son la solución.
#
# **Rango.** $\operatorname{rango}(A)$ es el número de pivotes de una forma escalonada de $A$, es decir, el número de renglones que no quedan en ceros al escalonar. Las columnas donde caen los pivotes son las **columnas pivote**.
#
# - No depende de qué operaciones elijas para escalonar.
# - Si $A$ es $m\times n$, $\operatorname{rango}(A)\le\min(m,n)$: cada pivote ocupa un renglón y una columna distintos.
# - Es el número máximo de renglones de $A$ entre los que ninguno es combinación de los otros, y coincide con el de columnas.
# - Si $A$ es $n\times n$: $\operatorname{rango}(A)=n$ si y solo si $\det A\neq 0$, si y solo si $A$ es invertible.
#
# **Dos rangos.** Con $r=\operatorname{rango}(A)$ y $r^*=\operatorname{rango}[A\,|\,\mathbf b]$, siempre $r\le r^*\le r+1$. Si $r^*=r$, $\mathbf b$ es combinación de las columnas y el sistema tiene solución. Si $r^*=r+1$, no lo es: al escalonar $[A\,|\,\mathbf b]$ aparece un renglón $0=c$ con $c\neq0$.
#
# `numpy.linalg.matrix_rank` calcula el rango en punto flotante, con una tolerancia; la calculadora lo compara con el rango exacto.

# %%
def _solo_gauss(pasos):
    """Los pasos de la fase de Gauss (hasta la forma escalonada), sin los de Jordan."""
    fases = [k for k, (tipo, _, _) in enumerate(pasos) if tipo == "fase"]
    return pasos[:fases[1]] if len(fases) > 1 else pasos


def rango_y_pivotes(A, pasos=False):
    """(forma escalonada exacta, columnas pivote numeradas desde 1, rango) de una matriz A de hasta 6×6."""
    A = matriz_limitada(A, "A")
    el = _elimina(A, A.cols, con_b=False)
    if pasos:
        muestra("A", A)
        _imprime_pasos(_solo_gauss(el["pasos"]))
    return el["E"], [c + 1 for c in el["pivotes"]], el["r"]


def combinacion(A, b):
    """¿Es b combinación lineal de las columnas de A? Diccionario con 'es' (True o False), los rangos r y r_aum,
    las formas escalonada y reducida de [A | b], los coeficientes (con parámetros s, t, ... si hay infinitas
    combinaciones), la particular y las direcciones, y, si no lo es, el renglón '0 = c' que lo prueba."""
    A, b = matriz_limitada(A, "A"), vector(b, "b")
    if b.rows != A.rows:
        raise ValueError(f"b tiene {_entradas(b.rows)} y A tiene {A.rows} renglones: deben coincidir")
    n = A.cols
    el = _elimina(A.row_join(b), n, con_b=True)
    res = {"es": el["imposible"] is None, "r": el["r"], "r_aum": el["r_aum"], "escalonada": el["E"],
           "reducida": el["R"], "pivotes": [c + 1 for c in el["pivotes"]], "coeficientes": None, "particular": None,
           "direcciones": [], "parametros": [], "renglon": None}
    if el["imposible"] is not None:
        f, c = el["imposible"]
        res["renglon"] = f"R{f} de la forma escalonada de [A | b] dice 0 = {c}"
    else:
        x, xp, dirs, ps = solucion_general(el["R"], el["pivotes"], n)
        res.update(coeficientes=x, particular=xp, direcciones=dirs, parametros=ps)
    return res


def _texto_combinacion(coef, ps=()):
    """'2·a₁ + 1·a₂' o, con parámetros, '(3 - 2s)·a₁ + s·a₂' (la constante primero)."""
    partes = []
    for j, c in enumerate(coef, 1):
        c = sp.sympify(c)
        texto = _paren(c) if c.is_number else (str(c) if c.is_Symbol else f"({_lineal(c, ps)})")
        partes.append(f"{texto}·a{_sub(j)}")
    return " + ".join(partes)


def _flecha(ax, origen, punta, color, estilo, ancho, texto, desplazamiento=(6, 6)):
    if np.hypot(*(np.asarray(punta, float) - np.asarray(origen, float))) > 1e-12:
        ax.annotate("", xy=punta, xytext=origen,
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=ancho, linestyle=estilo, shrinkA=0, shrinkB=0))
    ax.annotate(texto, xy=punta, xytext=desplazamiento, textcoords="offset points", fontsize=11)


def dibuja_columnas(A, b=None):
    """Lectura por columnas en el plano (A de 2×2): a₁, a₂, b y, si se puede, el camino x₁a₁ + x₂a₂ hasta b."""
    A = matriz_limitada(A, "A")
    if A.shape != (2, 2):
        raise ValueError(f"el dibujo de las columnas es para A de 2×2 (esta es {dims(A)})")
    a1, a2 = [_flotante(A[:, j]).ravel() for j in range(2)]
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    puntos = [np.zeros(2), a1, a2]
    _flecha(ax, (0, 0), a1, _COLORES[0], "-", 2, "a₁", (8, -14))
    _flecha(ax, (0, 0), a2, _COLORES[1], "--", 2, "a₂")
    titulo = "Columnas de A"
    if b is not None:
        res = combinacion(A, b)
        bf = _flotante(vector(b, "b")).ravel()
        puntos.append(bf)
        if res["es"]:
            c = res["particular"]                   # con infinitas combinaciones, la que tiene los parámetros en 0
            p1 = float(c[0]) * a1
            puntos.append(p1)
            ax.plot(*np.array([[0, 0], p1]).T, ":", color="0.35", lw=1.5)
            ax.plot(*np.array([p1, bf]).T, ":", color="0.35", lw=1.5, label="camino x₁a₁, luego + x₂a₂")
            ax.plot(*p1, "s", color="0.35", ms=6)
            titulo = f"b = {_texto_combinacion(c)}" + (" (una de infinitas)" if res["parametros"] else "")
        else:
            titulo = "b no es combinación de a₁ y a₂"
        if res["r"] == 1:                           # columnas en una misma recta: se dibuja esa recta
            d = a1 if np.hypot(*a1) > 0 else a2
            t = np.linspace(-1, 1, 2) * 2 * max(np.abs(np.vstack(puntos)).max(), 1) / np.hypot(*d)
            ax.plot(t * d[0], t * d[1], "-.", color="0.6", lw=1, label="recta de las combinaciones")
        _flecha(ax, (0, 0), bf, "black", "-", 3, "b")
    lim = 1.2 * max(1.0, float(np.abs(np.vstack(puntos)).max()))
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
    ax.axhline(0, color="black", lw=0.6); ax.axvline(0, color="black", lw=0.6)
    ax.grid(alpha=0.3); ax.set_xlabel("primera componente"); ax.set_ylabel("segunda componente")
    ax.set_title(titulo)
    if ax.get_legend_handles_labels()[0]:
        ax.legend(loc="lower right", fontsize=8)
    plt.show()

# %% [markdown]
# **Recorrido paso a paso.** Primero, $A\mathbf x$ como combinación de columnas. Después, el rango de $\begin{bmatrix}1&2&3\\2&4&7\\3&6&10\end{bmatrix}$: al escalonar, la columna 2 se queda sin pivote porque es el doble de la columna 1.

# %%
_A = matriz("2, 1; 1, 3")
_x = vector("2, 1")
muestra("A x", _A * _x)
muestra("2·a₁ + 1·a₂", 2 * _A[:, 0] + 1 * _A[:, 1])
print("Son el mismo vector: A x es la combinación de las columnas con los coeficientes de x.\n")
_E, _piv, _r = rango_y_pivotes("1, 2, 3; 2, 4, 7; 3, 6, 10", pasos=True)
muestra("forma escalonada", _E)
print(f"Columnas pivote: {_lista(_piv)} → rango = {_r}; numpy.linalg.matrix_rank da",
      np.linalg.matrix_rank(_flotante(matriz("1, 2, 3; 2, 4, 7; 3, 6, 10"))))
dibuja_columnas(_A, _A * _x)

# %% [markdown]
# **Calculadora.** Escribe $A$ por renglones y, si quieres, $\mathbf b$ (déjalo vacío para calcular solo el rango de $A$). Si $A$ es $2\times 2$, puedes dibujar la lectura por columnas.

# %%
def _informe_rango(A_txt, b_txt, n_cifras):
    """Imprime el análisis y devuelve (A, b) para el dibujo, o None si la entrada tiene un error."""
    try:
        A = matriz_limitada(A_txt, "A")
        b = vector(b_txt, "b") if b_txt.strip() else None
        res = combinacion(A, b) if b is not None else None
    except ValueError as err:
        print("Revisa la entrada:", err); return None
    muestra("A", A, n_cifras)
    E, piv, r = rango_y_pivotes(A)
    muestra("forma escalonada de A", E)
    r_np = int(np.linalg.matrix_rank(_flotante(A)))
    print(f"Columnas pivote: {_lista(piv) if piv else 'ninguna'} → rango(A) = {r}")
    print(f"numpy.linalg.matrix_rank(A) = {r_np}" + (" (coincide)" if r_np == r else
          " (NO coincide: numpy usa una tolerancia en punto flotante; el rango exacto es el de la forma escalonada)"))
    if b is None:
        return A, None
    Ab = A.row_join(b)
    muestra_aumentada("forma escalonada de [A | b]", res["escalonada"], A.cols)
    r_np2 = int(np.linalg.matrix_rank(_flotante(Ab)))
    print(f"rango[A | b] = {res['r_aum']}; numpy.linalg.matrix_rank([A | b]) = {r_np2}" +
          (" (coincide)" if r_np2 == res["r_aum"] else " (NO coincide: numpy usa una tolerancia en punto flotante)"))
    if res["es"]:
        print(f"rango[A | b] = rango(A) = {r}: b SÍ es combinación de las columnas de A.")
        if res["parametros"]:
            print(f"Hay infinitas combinaciones ({len(res['parametros'])} parámetro{'s' if len(res['parametros']) > 1 else ''}, "
                  f"{_coma(res['parametros'])}): b = {_texto_combinacion(res['coeficientes'], res['parametros'])}")
            print(f"Por ejemplo, con {'todos los parámetros en 0' if len(res['parametros']) > 1 else str(res['parametros'][0]) + ' = 0'}: "
                  f"b = {_texto_combinacion(res['particular'])}")
        else:
            print(f"Una sola combinación: b = {_texto_combinacion(res['coeficientes'])}")
        print("Comprobación: A·(coeficientes) = b" if A * res["particular"] == b else "Comprobación FALLÓ")
    else:
        print(f"rango[A | b] = {res['r_aum']} = rango(A) + 1: b NO es combinación de las columnas de A; "
              f"{res['renglon']}.")
    return A, b


def calculadora_rango(A_txt, b_txt, dibujar, n_cifras):
    with bloque():                 # el texto aparece completo, de una sola vez; el dibujo va después
        dibujo = _informe_rango(A_txt, b_txt, n_cifras)
    if dibujar and dibujo is not None and dibujo[0].shape == (2, 2):
        dibuja_columnas(*dibujo)


widgets.interact(calculadora_rango,
    A_txt=widgets.Textarea(value="1, 2\n2, 4", description="A =", continuous_update=False,
                           layout=widgets.Layout(width="340px", height="110px")),
    b_txt=widgets.Text(value="3, 6", description="b =", continuous_update=False),
    dibujar=widgets.Checkbox(value=True, description="dibujar columnas (A de 2×2)"),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False));

# %%
# Casos de prueba de la sección 3 (resultado conocido)
_A13 = de_columnas("1, 0, 1", "0, 1, 1")
PRUEBAS_3 = [
    ("1, 2, 3; 2, 4, 7; 3, 6, 10 → rango 2, columnas pivote [1, 3]",
     lambda: rango_y_pivotes("1, 2, 3; 2, 4, 7; 3, 6, 10")[1:] == ([1, 3], 2)),
    ("identidad 3×3 → rango 3", lambda: rango_y_pivotes(sp.eye(3))[2] == 3),
    ("matriz de ceros 2×3 → rango 0 (ninguna columna pivote)", lambda: rango_y_pivotes(sp.zeros(2, 3))[1:] == ([], 0)),
    ("1, 2; 2, 4; 3, 6 → rango 1", lambda: rango_y_pivotes("1, 2; 2, 4; 3, 6")[2] == 1),
    ("1, 0, 2, 1; 0, 1, 1, 3 → rango 2", lambda: rango_y_pivotes("1, 0, 2, 1; 0, 1, 1, 3")[2] == 2),
    ("a₁ = (1, 0, 1), a₂ = (0, 1, 1), b = (1, 1, 2) → b = 1·a₁ + 1·a₂",
     lambda: (lambda r: r["es"] and r["coeficientes"] == vector("1, 1") and not r["parametros"])(combinacion(_A13, "1, 1, 2"))),
    ("mismas columnas, b = (1, 1, 3) → no es combinación: rango A = 2, rango [A|b] = 3, renglón 0 = 1",
     lambda: (lambda r: not r["es"] and (r["r"], r["r_aum"]) == (2, 3) and "0 = 1" in r["renglon"])(combinacion(_A13, "1, 1, 3"))),
    ("0.1, 0.2; 0.3, 0.6 → rango 1 con fracciones exactas, igual que numpy.linalg.matrix_rank",
     lambda: (lambda E, piv, r: r == 1 and all(e.is_Rational for e in E)
                                and np.linalg.matrix_rank(_flotante(matriz("0.1, 0.2; 0.3, 0.6"))) == 1)(
         *rango_y_pivotes("0.1, 0.2; 0.3, 0.6"))),
]
corre_pruebas(PRUEBAS_3, 3)

# %% [markdown]
# **Contrasta (sección 3).** Calcula a mano el rango de $\begin{bmatrix}1&2&3\\0&1&1\\1&3&4\end{bmatrix}$ escalonando, y escribe el número.

# %%
mi_rango = None     # escribe un número entero, por ejemplo 3

if mi_rango is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    try:
        mio = numero(mi_rango, "tu rango")
    except ValueError as err:
        print("Revisa tu número:", err)
    else:
        if not (mio.is_integer and 0 <= mio <= 3):     # respuesta imposible: la referencia todavía no se muestra
            print(f"El rango de una matriz 3×3 es un entero de 0 a 3 (escribiste {mio}): revisa tu número "
                  "y vuelve a ejecutar la celda.")
        else:
            E, piv, ref = rango_y_pivotes("1, 2, 3; 0, 1, 1; 1, 3, 4")
            print("La calculadora da: rango =", ref)
            print("coinciden" if mio == ref else
                  "NO coinciden: cuenta los renglones que no quedan en ceros después de escalonar; "
                  "fíjate si un renglón es combinación de los otros.")
            print("\nEl procedimiento completo:")
            rango_y_pivotes("1, 2, 3; 0, 1, 1; 1, 3, 4", pasos=True)
            muestra("forma escalonada", E)
