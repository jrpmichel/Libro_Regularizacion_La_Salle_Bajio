# ID: ALG-U3-NB02
# Notebook: alg/u3_sistemas.ipynb · sección 3.2 regla de Cramer y su costo
# Repositorio: alg/u3_sistemas/02_cramer.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Regla de Cramer
#
# Si $A$ es $n\times n$ y $\det A\neq 0$, el sistema $A\mathbf x=\mathbf b$ tiene una sola solución, y cada incógnita es un cociente de determinantes:
#
# $$x_i=\frac{\det A_i}{\det A},\qquad i=1,\dots,n,$$
#
# donde $A_i$ es la matriz $A$ con su **columna $i$ cambiada por $\mathbf b$**.
#
# **Por qué funciona.** Si $\mathbf b=x_1\mathbf a_1+\dots+x_n\mathbf a_n$ (con $\mathbf a_j$ la columna $j$ de $A$), en $A_i$ puedes restar a la columna $i$ los múltiplos $x_j\mathbf a_j$ de las demás columnas sin cambiar el determinante. Queda $x_i\mathbf a_i$ en su lugar y, al sacar el factor $x_i$, $\det A_i=x_i\det A$.
#
# **Si $\det A=0$, no se divide.** La regla no da solución; solo asegura que no hay solución única. Hay dos casos:
#
# - Algún $\det A_i\neq 0$: **no hay solución**, porque cualquier solución obligaría a $\det A_i=x_i\cdot\det A=x_i\cdot 0=0$.
# - Todos los $\det A_i$ valen 0: **Cramer no decide**. El sistema puede tener infinitas soluciones o ninguna; por ejemplo, tres planos paralelos distintos tienen todos los determinantes en cero y no se cortan. Se decide comparando rangos (sección 4).
#
# **Costo.** Con cofactores y sin aprovechar ceros, un determinante $n\times n$ cuesta $M_n=n\,(M_{n-1}+1)$ multiplicaciones, con $M_2=2$. Cramer calcula $n+1$ determinantes y hace $n$ divisiones: $(n+1)M_n+n$ operaciones. La eliminación de Gauss, contando multiplicaciones y divisiones, hace $\frac{n^3+3n^2-n}{3}$:
#
# | $n$ | Cramer con cofactores | Gauss |
# |---:|---:|---:|
# | 2 | 8 | 6 |
# | 3 | 39 | 17 |
# | 4 | 204 | 36 |
# | 10 | 68 588 310 | 430 |
#
# Cramer crece como $n!$ y Gauss como $n^3$. Úsalo en sistemas $2\times 2$ y $3\times 3$, o cuando quieras una fórmula con letras para una incógnita; para sistemas más grandes, usa eliminación.

# %%
def cramer(Ab, pasos=True, n_cifras=4):
    """Regla de Cramer exacta para un sistema cuadrado [A | b].

    Devuelve un diccionario con det A, la lista de det Aᵢ, las matrices Aᵢ, la solución (o None), el veredicto
    ('única', 'sin solución' o 'usar el rango') y un mensaje. Si det A = 0 no divide."""
    Ab = aumentada(Ab)
    m, n = Ab.rows, Ab.cols - 1
    if m != n:
        raise ValueError(f"la regla de Cramer solo se aplica a sistemas cuadrados (tantas ecuaciones como incógnitas); "
                         f"este tiene {_ecuaciones(m)} y {_incognitas(n)}: resuélvelo con Gauss-Jordan (sección 1) "
                         "o clasifícalo con los rangos (sección 4)")
    A, b = Ab[:, :n], Ab[:, n]
    d = A.det()
    matrices, dets = [], []
    for i in range(n):
        Ai = A.copy()
        Ai[:, i] = b
        matrices.append(Ai)
        dets.append(Ai.det())
    nombres = nombres_incognitas(n)
    if d != 0:
        x = sp.Matrix([di / d for di in dets])
        veredicto, mensaje = "única", f"det A = {d} ≠ 0: el sistema tiene una sola solución"
    else:
        x = None
        distintos = [i for i in range(n) if dets[i] != 0]
        if distintos:
            i = distintos[0]
            veredicto = "sin solución"
            mensaje = (f"det A = 0 y det A{_sub(i + 1)} = {dets[i]} ≠ 0: no hay solución. Si la hubiera, "
                       f"det A{_sub(i + 1)} = {nombres[i]}·det A = {nombres[i]}·0 = 0, y no es 0")
        else:
            veredicto = "usar el rango"
            mensaje = ("det A = 0 y todos los det Aᵢ valen 0: con Cramer no se puede decidir. El sistema puede tener "
                       "infinitas soluciones o ninguna; hay que usar el rango: compara rango(A) con rango[A | b] (sección 4)")
    res = {"det_A": d, "det_Ai": dets, "Ai": matrices, "solucion": x, "veredicto": veredicto, "mensaje": mensaje}
    if pasos:
        _informe_cramer(Ab, res, n_cifras)
    return res


def _informe_cramer(Ab, res, n_cifras):
    n = Ab.cols - 1
    A, b = Ab[:, :n], Ab[:, n]
    nombres = nombres_incognitas(n)
    _imprime_sistema(Ab)
    muestra("A", A)
    muestra("b", b)
    muestra("det A", res["det_A"], n_cifras)
    for i, (Ai, di) in enumerate(zip(res["Ai"], res["det_Ai"]), 1):
        print(f"A{_sub(i)}: la columna {i} de A cambiada por b")
        muestra(f"A{_sub(i)}", Ai)
        muestra(f"det A{_sub(i)}", di, n_cifras)
    print(res["mensaje"] + ".")
    if res["solucion"] is None:
        return
    d = res["det_A"]
    for i, (v, di, xi) in enumerate(zip(nombres, res["det_Ai"], res["solucion"]), 1):
        decimal = "" if xi.is_integer else f" ≈ {cifras(xi, n_cifras)}"
        print(f"  {v} = det A{_sub(i)} / det A = {_paren(di)}/{_paren(d)} = {xi}{decimal}")
    print("Comprobación: A x = b" if A * res["solucion"] == b else "Comprobación FALLÓ: A x ≠ b")


def costos(n):
    """(Cramer con cofactores, Gauss): multiplicaciones y divisiones para un sistema n×n, sin aprovechar ceros.
    Cramer: (n+1)·Mₙ + n con Mₙ = n·(Mₙ₋₁ + 1), M₂ = 2 (y M₁ = 0). Gauss: (n³ + 3n² - n)/3, que siempre es entero."""
    if isinstance(n, bool) or not isinstance(n, (int, np.integer, sp.Integer)) or int(n) < 1:
        raise ValueError(f"n = {n} no sirve: el tamaño del sistema debe ser un entero positivo")
    n = int(n)
    M = 0                                   # M₁ = 0: el determinante de una 1×1 es su única entrada
    for k in range(2, n + 1):
        M = k * (M + 1)                     # M₂ = 2, M₃ = 9, M₄ = 40, ...
    gauss, resto = divmod(n ** 3 + 3 * n ** 2 - n, 3)
    assert resto == 0
    return (n + 1) * M + n, gauss


def tabla_costos(n_sistema=None, tamanos=(2, 3, 4, 5, 6, 10)):
    """Imprime la tabla de operaciones de Cramer y de Gauss; marca el tamaño de tu sistema."""
    filas = sorted(set(tamanos) | ({n_sistema} if n_sistema else set()))
    print("Operaciones (multiplicaciones y divisiones, sin aprovechar ceros):")
    print(f"{'n':>4}  {'Cramer con cofactores':>22}  {'Gauss':>7}")
    for k in filas:
        c, g = costos(k)
        marca = "   ← tu sistema" if k == n_sistema else ""
        print(f"{k:>4}  {c:>22,}  {g:>7,}{marca}".replace(",", " "))

# %% [markdown]
# **Recorrido paso a paso.** El sistema $2x+y=5$, $x+3y=5$: primero $\det A$, luego cada $A_i$ (la columna $i$ cambiada por $\mathbf b$) y su determinante, y al final los cocientes. Después, la tabla de costos.

# %%
_recorrido = cramer("2, 1, 5; 1, 3, 5")
print()
tabla_costos()

# %% [markdown]
# **Calculadora.** Escribe un sistema cuadrado como matriz aumentada (la última columna es $\mathbf b$). Prueba también uno con $\det A=0$, por ejemplo `1, 2, 3; 2, 4, 7`.

# %%
def calculadora_cramer(Ab_txt, n_cifras):
    with bloque():                 # el resultado aparece completo, de una sola vez
        try:
            Ab = aumentada(Ab_txt)
            res = cramer(Ab, pasos=False)
        except ValueError as err:
            print("Revisa la entrada:", err); return
        _informe_cramer(Ab, res, n_cifras)
        print()
        tabla_costos(Ab.rows)


widgets.interact(calculadora_cramer,
    Ab_txt=widgets.Textarea(value="3, 2, 7\n1, -1, -1", description="[A | b] =", continuous_update=False,
                            layout=widgets.Layout(width="340px", height="110px")),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False));

# %%
# Casos de prueba de la sección 2 (resultado conocido)
_CR = lambda texto: cramer(texto, pasos=False)
PRUEBAS_2 = [
    ("2, 1, 5; 1, 3, 5 → det A = 5, det A₁ = 10, det A₂ = 5, x = (2, 1)",
     lambda: (lambda r: r["det_A"] == 5 and r["det_Ai"] == [10, 5] and r["solucion"] == vector("2, 1"))(_CR("2, 1, 5; 1, 3, 5"))),
    ("3, 2, 7; 1, -1, -1 → det A = -5, x = (1, 2)",
     lambda: (lambda r: r["det_A"] == -5 and r["solucion"] == vector("1, 2"))(_CR("3, 2, 7; 1, -1, -1"))),
    ("1, 1, 1, 6; 2, 2, 3, 15; 1, 3, 2, 13 → det A = -2, det A₃ = -6, x = (1, 2, 3)",
     lambda: (lambda r: r["det_A"] == -2 and r["det_Ai"][2] == -6 and r["solucion"] == vector("1, 2, 3"))(
         _CR("1, 1, 1, 6; 2, 2, 3, 15; 1, 3, 2, 13"))),
    ("1, 2, 3; 2, 4, 7 → det A = 0 y det A₁ = -2: no hay solución (sin dividir entre 0)",
     lambda: (lambda r: r["det_A"] == 0 and r["det_Ai"][0] == -2 and r["solucion"] is None
                        and r["veredicto"] == "sin solución" and "no hay solución" in r["mensaje"])(_CR("1, 2, 3; 2, 4, 7"))),
    ("1, 2, 3; 2, 4, 6 → det A = 0 y todos los det Aᵢ = 0: hay que usar el rango",
     lambda: (lambda r: r["det_A"] == 0 and r["det_Ai"] == [0, 0] and r["veredicto"] == "usar el rango"
                        and "hay que usar el rango" in r["mensaje"])(_CR("1, 2, 3; 2, 4, 6"))),
    ("planos paralelos 1, 1, 1, 1; 1, 1, 1, 2; 1, 1, 1, 3 → todos los determinantes 0: hay que usar el rango",
     lambda: (lambda r: r["det_A"] == 0 and r["det_Ai"] == [0, 0, 0] and r["veredicto"] == "usar el rango"
                        and "hay que usar el rango" in r["mensaje"])(_CR("1, 1, 1, 1; 1, 1, 1, 2; 1, 1, 1, 3"))),
    ("costos(2) = (8, 6), costos(3) = (39, 17), costos(4) = (204, 36), costos(10) = (68588310, 430)",
     lambda: [costos(k) for k in (2, 3, 4, 10)] == [(8, 6), (39, 17), (204, 36), (68588310, 430)]),
    ("0.5, 0.3, 1.9; 0.2, 0.4, 1.4 → x = (17/7, 16/7) exacto",
     lambda: _CR("0.5, 0.3, 1.9; 0.2, 0.4, 1.4")["solucion"] == sp.Matrix([sp.Rational(17, 7), sp.Rational(16, 7)])),
]
corre_pruebas(PRUEBAS_2, 2)

# %% [markdown]
# **Contrasta (sección 2).** Resuelve con Cramer $2x+5y=1$, $x+3y=2$: calcula $\det A$, $\det A_1$ y $\det A_2$ y divide. Escribe tu solución (primero $x$, luego $y$) y, si quieres, tu $\det A$.

# %%
mi_solucion = None     # escribe tu solución como texto, primero x y luego y; por ejemplo: "4, -1"
mi_det_A = None        # opcional: tu det A, por ejemplo "-3"

if mi_solucion is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    try:
        mia = vector(str(mi_solucion), "tu solución")
        mio_det = None if mi_det_A is None else numero(mi_det_A, "tu det A")
    except ValueError as err:
        print("Revisa tu respuesta:", err)
    else:
        if mia.rows != 2:          # respuesta incompleta: la referencia todavía no se muestra
            print(f"Tu solución debe tener dos números, primero x y luego y (escribiste {mia.rows}): "
                  "complétala y vuelve a ejecutar la celda.")
        else:
            ref = cramer("2, 5, 1; 1, 3, 2", pasos=False)
            print(f"La calculadora da: det A = {ref['det_A']}; x = {ref['solucion'][0]}, y = {ref['solucion'][1]}")
            if mio_det is not None:
                print("det A:", "coincide" if mio_det == ref["det_A"] else
                      "NO coincide: det A = a₁₁a₂₂ - a₁₂a₂₁ (el producto de la diagonal menos el de la otra diagonal)")
            print("Solución:", compara_vector(mia, ref["solucion"], ["x", "y"],
                                              "Para det A₁ cambia la columna de x por b = (1, 2); para det A₂, la de y."))
            print("\nEl procedimiento completo:")
            cramer("2, 5, 1; 1, 3, 2")
