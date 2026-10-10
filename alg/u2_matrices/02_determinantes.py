# ID: ALG-U2-NB02
# Notebook: alg/u2_matrices.ipynb · sección 2.2 determinantes por Sarrus y por cofactores
# Repositorio: alg/u2_matrices/02_determinantes.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Determinantes: Sarrus y cofactores
#
# El determinante es un número asociado a una matriz **cuadrada**. Para una $2\times 2$:
#
# $$\det\begin{bmatrix}a&b\\c&d\end{bmatrix}=ad-bc.$$
#
# **Menores y cofactores.** El menor $M_{ij}$ es el determinante de la submatriz que queda al borrar el renglón $i$ y la columna $j$. El cofactor le agrega un signo que alterna como un tablero de ajedrez:
#
# $$C_{ij}=(-1)^{i+j}M_{ij},\qquad \text{signos: }\begin{bmatrix}+&-&+\\-&+&-\\+&-&+\end{bmatrix}.$$
#
# **Desarrollo por cofactores.** Se hace por cualquier renglón $i$ o por cualquier columna $j$:
#
# $$\det A=\sum_{j=1}^{n}a_{ij}C_{ij}\quad(\text{renglón } i)\qquad\text{o}\qquad \det A=\sum_{i=1}^{n}a_{ij}C_{ij}\quad(\text{columna } j).$$
#
# Todas las líneas dan el mismo número, así que conviene la que tiene más ceros: cada $a_{ij}=0$ elimina un término completo, con todo y su cofactor.
#
# **Regla de Sarrus (solo $3\times 3$).** Copia las dos primeras columnas a la derecha de la matriz; suma los productos de las tres diagonales que bajan y resta los de las tres que suben:
#
# $$\det A=a_{11}a_{22}a_{33}+a_{12}a_{23}a_{31}+a_{13}a_{21}a_{32}-a_{13}a_{22}a_{31}-a_{11}a_{23}a_{32}-a_{12}a_{21}a_{33}.$$
#
# No se generaliza: un determinante $4\times 4$ tiene $4!=24$ términos y las "diagonales" de Sarrus solo darían 8.
#
# **Costo.** Sin ceros, el desarrollo por cofactores de una $n\times n$ hace $\sum_{k=1}^{n-1} n!/k!\approx 1.7\,n!$ multiplicaciones: 9 para $n=3$, 40 para $n=4$, 205 para $n=5$ y unos 6.2 millones para $n=10$. Por eso el software no calcula determinantes así, sino por eliminación (sección 5), que hace del orden de $n^3/3$.
#
# La calculadora muestra qué renglón o columna usa en cada nivel, cada término $a_{ij}C_{ij}$, cuántas multiplicaciones hizo (un término con $a_{ij}=0$ no cuenta) y lo compara con el determinante exacto de `sympy`.

# %%
def _indice(k, n, letra):
    if isinstance(k, bool) or not isinstance(k, (int, np.integer, sp.Integer)) or not 1 <= int(k) <= n:
        raise ValueError(f"{letra} = {k} no sirve: i y j van de 1 a {n} (la matriz es {n}×{n})")
    return int(k)


def submatriz(A, i, j):
    """A sin el renglón i ni la columna j (i, j desde 1)."""
    A = como_matriz(A, "A")
    _cuadrada(A, "los menores solo existen para matrices cuadradas")
    if A.rows < 2:
        raise ValueError("una matriz 1×1 no tiene menores: su determinante es su única entrada")
    i, j = _indice(i, A.rows, "i"), _indice(j, A.rows, "j")
    return A.minor_submatrix(i - 1, j - 1)


def menor(A, i, j):
    """M_ij: determinante de la submatriz sin el renglón i ni la columna j."""
    return submatriz(A, i, j).det()


def cofactor(A, i, j):
    """C_ij = (-1)^(i+j) M_ij."""
    return (-1) ** (int(i) + int(j)) * menor(A, i, j)


_BAJAN = [((1, 1), (2, 2), (3, 3)), ((1, 2), (2, 3), (3, 1)), ((1, 3), (2, 1), (3, 2))]
_SUBEN = [((1, 3), (2, 2), (3, 1)), ((1, 1), (2, 3), (3, 2)), ((1, 2), (2, 1), (3, 3))]


def det_sarrus(A):
    """(det A, texto con los pasos) por la regla de Sarrus; solo para 3×3."""
    A = como_matriz(A, "A")
    if A.shape != (3, 3):
        raise ValueError(f"la regla de Sarrus solo vale para matrices 3×3 (esta es {dims(A)}); usa cofactores")
    a = lambda p: A[p[0] - 1, p[1] - 1]
    cuenta = lambda d: "·".join(_paren(a(p)) for p in d)
    bajan = [a(d[0]) * a(d[1]) * a(d[2]) for d in _BAJAN]
    suben = [a(d[0]) * a(d[1]) * a(d[2]) for d in _SUBEN]
    det = sum(bajan) - sum(suben)
    pasos = "\n".join([
        "Copia las dos primeras columnas a la derecha de A:",
        _plano(A.row_join(A[:, :2]), corte=3),
        f"Diagonales que bajan (+): {' + '.join(cuenta(d) for d in _BAJAN)} = {_suma(bajan)} = {sum(bajan)}",
        f"Diagonales que suben (-): {' + '.join(cuenta(d) for d in _SUBEN)} = {_suma(suben)} = {sum(suben)}",
        f"det A = {sum(bajan)} - {_paren(sum(suben))} = {det}",
    ])
    return det, pasos


def _mejor_linea(A):
    """('renglón' o 'columna', índice desde 1, ceros) de la línea con más ceros; en empate gana la primera."""
    mejor = ("renglón", 1, -1)
    for tipo, lineas in (("renglón", [A.row(i) for i in range(A.rows)]), ("columna", [A.col(j) for j in range(A.cols)])):
        for k, linea in enumerate(lineas, 1):
            ceros = sum(1 for e in linea if e == 0)
            if ceros > mejor[2]:
                mejor = (tipo, k, ceros)
    return mejor


def _expande(A, nivel, pasos, nombre):
    """(det, multiplicaciones) por cofactores; si pasos es una lista, le agrega el texto de cada paso."""
    n, s = A.rows, "    " * nivel
    if n == 1:
        if pasos is not None and nivel == 0:
            pasos.append(f"det {nombre} = {A[0, 0]} (una matriz 1×1: su determinante es su única entrada)")
        return A[0, 0], 0
    if n == 2:
        a, b, c, d = A
        det = a * d - b * c
        if pasos is not None:
            pasos.append(f"{s}det {nombre} = {_paren(a)}·{_paren(d)} - {_paren(b)}·{_paren(c)} = {det}")
        return det, int(a != 0 and d != 0) + int(b != 0 and c != 0)
    tipo, k, ceros = _mejor_linea(A)
    if pasos is not None:
        nota = f" (tiene {ceros} cero{'s' if ceros != 1 else ''})" if ceros else ""
        indices = " (índices de la submatriz)" if nivel > 0 else ""      # a_ij y C_ij de abajo cuentan dentro del menor
        pasos.append(f"{s}det {nombre} ({dims(A)}): desarrollo por el {tipo} {k}{nota}{indices}")
    total, mult, terminos = 0, 0, []
    for t in range(1, n + 1):
        i, j = (k, t) if tipo == "renglón" else (t, k)
        a = A[i - 1, j - 1]
        if a == 0:
            if pasos is not None:
                pasos.append(f"{s}  a{_sub(i, j)} = 0: ese término vale 0 y su cofactor no se calcula")
            continue
        signo = "+1" if (i + j) % 2 == 0 else "-1"
        sub = A.minor_submatrix(i - 1, j - 1)
        if pasos is not None:
            pasos.append(f"{s}  a{_sub(i, j)}·C{_sub(i, j)} = {_paren(a)}·({signo})·M{_sub(i, j)},  M{_sub(i, j)} = det {_en_linea(sub)}")
        M, m = _expande(sub, nivel + 1, pasos, _en_linea(sub))
        termino = a * (-1) ** (i + j) * M
        mult += m + 1                                  # las del menor y la de a_ij·C_ij
        total += termino
        terminos.append(termino)
        if pasos is not None:
            pasos.append(f"{s}  a{_sub(i, j)}·C{_sub(i, j)} = {_paren(a)}·({signo})·{_paren(M)} = {termino}")
    if pasos is not None:
        if len(terminos) > 1:
            pasos.append(f"{s}det {nombre} = {_suma(terminos)} = {total}")
        elif terminos:
            pasos.append(f"{s}det {nombre} = {total} (solo un término no es cero)")
        else:
            pasos.append(f"{s}det {nombre} = 0: todas las entradas de esa línea son cero")
    return total, mult


def det_cofactores(A, pasos=False):
    """(det A, lista de pasos, multiplicaciones) por desarrollo de cofactores, siempre por la línea con más ceros."""
    A = como_matriz(A, "A")
    _cuadrada(A)
    lista = [] if pasos else None
    det, mult = _expande(A, 0, lista, "A")
    return det, (lista or []), mult


def multiplicaciones_sin_ceros(n):
    """Multiplicaciones del desarrollo por cofactores de una n×n sin ceros: f(n) = n(1 + f(n-1)), f(2) = 2."""
    f = 0
    for k in range(2, n + 1):
        f = k * (1 + f) if k > 2 else 2
    return f


_MAX_COFACTORES = 7      # 8! = 40320 términos: más que eso ya tarda demasiado
_MAX_PASOS = 80          # renglones de pasos que se imprimen completos


def calculadora_determinante(A_txt, metodo, n_cifras):
    try:
        A = matriz(A_txt, "A")
        _cuadrada(A)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    muestra("A", A, n_cifras)
    exacto = A.det()
    calculados = []
    if metodo in ("Sarrus", "los dos"):
        try:
            d, texto = det_sarrus(A)
        except ValueError as err:
            if metodo == "Sarrus":
                print("Revisa la entrada:", err); return
            print("Sarrus no aplica:", err)
        else:
            print("— Regla de Sarrus —")
            print(texto)
            print("Multiplicaciones: 12 (6 productos de 3 factores, aunque haya ceros)")
            calculados.append(("Sarrus", d))
    if metodo in ("cofactores", "los dos"):
        if A.rows > _MAX_COFACTORES:
            print(f"Cofactores: con una {dims(A)} serían hasta {multiplicaciones_sin_ceros(A.rows)} multiplicaciones; "
                  f"esta calculadora los desarrolla hasta {_MAX_COFACTORES}×{_MAX_COFACTORES}. Usa la factorización LU (sección 5).")
            print(f"det A (sympy) = {exacto}" + ("" if exacto.is_integer else f" ≈ {cifras(exacto, n_cifras)}"))
        else:
            d, pasos, mult = det_cofactores(A, pasos=True)
            print("— Desarrollo por cofactores —")
            if len(pasos) <= _MAX_PASOS:
                print("\n".join(pasos))
            else:
                print("\n".join(p for p in pasos if not p.startswith("    ")))
                print(f"(se omiten {sum(1 for p in pasos if p.startswith('    '))} renglones de los menores internos)")
            print(f"Multiplicaciones: {mult} (una {dims(A)} sin ceros pediría {multiplicaciones_sin_ceros(A.rows)})")
            calculados.append(("cofactores", d))
    for metodo_usado, d in calculados:
        print(f"{metodo_usado}: det A = {d}" + ("" if d.is_integer else f" ≈ {cifras(d, n_cifras)}") +
              (" | coincide con sympy" if d == exacto else f" | NO coincide con sympy ({exacto})"))


widgets.interact(calculadora_determinante,
    A_txt=widgets.Textarea(value="2, -1, 3\n0, 4, 1\n5, 2, -2", description="A =", continuous_update=False,
                           layout=widgets.Layout(width="320px", height="110px")),
    metodo=widgets.Dropdown(options=["cofactores", "Sarrus", "los dos"], value="los dos", description="método"),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False));

# %%
# Casos de prueba de la sección 2 (resultado conocido)
_D3 = matriz("2, -1, 3; 0, 4, 1; 5, 2, -2")
_D4 = matriz("1, 0, 2, 0; 0, 3, 0, 1; 2, 0, 1, 0; 0, 1, 0, 2")
PRUEBAS_2 = [
    ("det [[3,1],[4,2]] = 2", lambda: det_cofactores(matriz("3, 1; 4, 2"))[0] == 2),
    ("det [[2,-1,3],[0,4,1],[5,2,-2]] = -85 por Sarrus y por cofactores",
     lambda: det_sarrus(_D3)[0] == -85 and det_cofactores(_D3)[0] == -85),
    ("en esa matriz, M₁₂ = -5 y C₁₂ = 5", lambda: menor(_D3, 1, 2) == -5 and cofactor(_D3, 1, 2) == 5),
    ("4×4 por cofactores = det de sympy (-15)", lambda: det_cofactores(_D4)[0] == _D4.det() == -15),
    ("3×3 sin ceros: 9 multiplicaciones; con un cero, menos",
     lambda: det_cofactores(matriz("1, 2, 3; 4, 5, 6; 7, 8, 10"))[2] == 9 == multiplicaciones_sin_ceros(3)
             and det_cofactores(_D3)[2] < 9),
    ("Sarrus en una 4×4 se rechaza", lambda: "solo vale para matrices 3×3" in _mensaje(lambda: det_sarrus(_D4))),
    ("una matriz no cuadrada se rechaza", lambda: _rechaza(lambda: det_cofactores(matriz("1, 2, 3; 4, 5, 6")))),
    ("i o j fuera de rango se rechazan", lambda: _rechaza(lambda: menor(_D3, 4, 1)) and _rechaza(lambda: cofactor(_D3, 0, 2))),
]
for nombre, prueba in PRUEBAS_2:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 2).** Calcula a mano $\det\begin{bmatrix}1&2&0\\3&-1&2\\0&1&1\end{bmatrix}$ por Sarrus y por cofactores; los dos deben darte el mismo número. Escríbelo.

# %%
mi_det = None     # escribe un número (también sirve como texto, por ejemplo "-7")

if mi_det is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    try:
        mio = numero(mi_det, "tu determinante")
    except ValueError as err:
        print("Revisa tu número:", err)
    else:
        ref = det_cofactores(matriz("1, 2, 0; 3, -1, 2; 0, 1, 1"))[0]
        print("La calculadora da: det =", ref)
        print("coinciden" if mio == ref else
              "NO coinciden: desarrolla por el renglón 1 o la columna 1 (tienen un cero) y cuida el signo (-1)^(i+j) de cada cofactor")
