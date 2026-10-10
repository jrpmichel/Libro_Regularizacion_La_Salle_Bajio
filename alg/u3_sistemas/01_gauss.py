# ID: ALG-U3-NB01
# Notebook: alg/u3_sistemas.ipynb · sección 3.1 eliminación de Gauss y Gauss-Jordan
# Repositorio: alg/u3_sistemas/01_gauss.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Eliminación de Gauss y Gauss-Jordan
#
# Un sistema de $m$ ecuaciones con $n$ incógnitas se guarda en su **matriz aumentada** $[A\,|\,\mathbf b]$: cada renglón es una ecuación, las primeras $n$ columnas son los coeficientes y la última es $\mathbf b$. Por ejemplo,
#
# $$\begin{aligned}2x+3y&=8\\ x-y&=-1\end{aligned}\qquad\longleftrightarrow\qquad\left[\begin{array}{cc|c}2&3&8\\1&-1&-1\end{array}\right].$$
#
# **Operaciones elementales de renglón.** No cambian las soluciones, porque cada una se deshace con otra del mismo tipo:
#
# - intercambiar dos renglones: $R_i\leftrightarrow R_j$;
# - multiplicar un renglón por un número $k\neq 0$: $R_i\leftarrow k\,R_i$;
# - sumar a un renglón un múltiplo de **otro**: $R_i\leftarrow R_i-k\,R_j$, con $i\neq j$.
#
# Cada operación actúa sobre el renglón completo, **incluida la columna de $\mathbf b$**.
#
# **Formas escalonadas.** Una matriz está en *forma escalonada* si los renglones de ceros van al final y el **pivote** de cada renglón (su primera entrada distinta de cero) queda a la derecha del pivote del renglón de arriba. Está en *forma escalonada reducida* si, además, cada pivote vale 1 y es lo único distinto de cero en su columna.
#
# - **Gauss** baja hasta la forma escalonada: en la columna del pivote $a_{kk}$, a cada renglón $i$ de abajo le restas $m_{ik}=a_{ik}/a_{kk}$ veces el renglón $k$, $R_i\leftarrow R_i-m_{ik}R_k$. Después despejas de la última ecuación a la primera (sustitución hacia atrás).
# - **Gauss-Jordan** sigue hasta la forma reducida: cada pivote en 1 y ceros también arriba de él. Si la solución es única, queda escrita en la última columna.
#
# Tres situaciones cambian el camino:
#
# 1. **Pivote cero.** Si en la posición del pivote hay un 0 y abajo hay una entrada distinta de cero, intercambias esos renglones.
# 2. **Columna sin pivote.** Si de ese renglón hacia abajo toda la columna es cero, esa columna no tiene pivote y pasas a la siguiente. El sistema ya no tiene solución única; la sección 4 dice si tiene infinitas o ninguna.
# 3. **Renglón $(0\ \cdots\ 0\,|\,c)$ con $c\neq 0$.** Dice $0=c$, que nadie cumple: el sistema no tiene solución.
#
# **Inversa por Gauss-Jordan.** Si $A$ es $n\times n$, Gauss-Jordan aplicado a $[A\,|\,I_n]$ termina en $[I_n\,|\,A^{-1}]$. Si alguna columna de $A$ se queda sin pivote, $A$ no tiene inversa (su determinante es 0).
#
# **Comprueba siempre en las ecuaciones originales.** Si un renglón se operó mal, la solución cumple la ecuación equivocada, y solo la sustitución en el sistema de partida lo detecta: el **residuo** $\mathbf b-A\mathbf x$ debe ser $\mathbf 0$.
#
# La calculadora trabaja con fracciones exactas e intercambia renglones solo cuando el pivote es cero, como harías a mano. Al final compara con `numpy.linalg.solve`, que trabaja en punto flotante y además intercambia renglones para usar el pivote más grande.

# %%
def _filas_independientes(A):
    """Índices (desde 0) de renglones de A que no son combinación de los anteriores."""
    elegidas = []
    for i in range(A.rows):
        if sp.Matrix.vstack(*[A.row(k) for k in elegidas + [i]]).rank() > len(elegidas):
            elegidas.append(i)
    return elegidas


def gauss_jordan(Ab, pasos=True, n_cifras=4):
    """Resuelve el sistema [A | b] por Gauss-Jordan con fracciones exactas.

    Devuelve un diccionario con la forma escalonada, la reducida, las operaciones de renglón (texto), las columnas con
    y sin pivote (desde 1), si la solución es única, la solución, el residuo b - A x y la solución de numpy.linalg.solve
    (solo si la solución es única; si numpy no puede, None y un aviso). Con pasos=True imprime cada operación."""
    Ab = aumentada(Ab)
    m, n = Ab.rows, Ab.cols - 1
    A, b = Ab[:, :n], Ab[:, n]
    nombres = nombres_incognitas(n)
    el = _elimina(Ab, n, con_b=True)
    unica = el["r"] == n and el["imposible"] is None
    res = {"escalonada": el["E"], "reducida": el["R"], "operaciones": operaciones(el), "pasos": el["pasos"],
           "pivotes": [c + 1 for c in el["pivotes"]], "sin_pivote": [c + 1 for c in el["libres"]],
           "imposible": el["imposible"], "unica": unica, "solucion": None, "residuo": None,
           "numpy": None, "aviso_numpy": None, "ecuaciones_numpy": None, "mensajes": []}
    if el["libres"]:
        cols = el["libres"]
        cuales = "la columna" if len(cols) == 1 else "las columnas"
        texto = f"{cuales} {_lista([c + 1 for c in cols])} no {'tiene' if len(cols) == 1 else 'tienen'} pivote"
        if el["imposible"] is None:
            libres = _coma([nombres[c] for c in cols])
            texto += (f" ({libres} {'es libre' if len(cols) == 1 else 'son libres'}): no hay solución única; "
                      "la clasificación (infinitas soluciones o ninguna) se hace en la sección 4")
        res["mensajes"].append(texto)
    if el["imposible"] is not None:
        f, c = el["imposible"]
        res["mensajes"].append(f"aparece el renglón 0 = {c} (R{f} de la forma escalonada): el sistema no tiene solución; "
                               "la clasificación se hace en la sección 4")
    if unica:
        x = el["R"][:n, n]
        res["solucion"] = x
        res["residuo"] = b - A * x
        filas = list(range(m)) if m == n else _filas_independientes(A)
        res["numpy"], res["aviso_numpy"] = _numpy_solve(A.extract(filas, list(range(n))), b.extract(filas, [0]))
        res["ecuaciones_numpy"] = [i + 1 for i in filas]
    if pasos:
        _informe_gauss(Ab, res, n_cifras, con_operaciones=True)
    return res


def _sustitucion(fila, x):
    """'2·1 + 3·2' para el renglón (2, 3) y x = (1, 2)."""
    return " + ".join(f"{_paren(a)}·{_paren(v)}" for a, v in zip(fila, x))


def _informe_gauss(Ab, res, n_cifras, con_operaciones=True):
    m, n = Ab.rows, Ab.cols - 1
    A, b = Ab[:, :n], Ab[:, n]
    nombres = nombres_incognitas(n)
    _imprime_sistema(Ab)
    muestra_aumentada("[A | b]", Ab, n)
    if con_operaciones:
        _imprime_pasos(res["pasos"], corte=n)
    muestra_aumentada("forma escalonada", res["escalonada"], n)
    muestra_aumentada("forma escalonada reducida", res["reducida"], n)
    print("Columnas con pivote:", _lista(res["pivotes"]) if res["pivotes"] else "ninguna")
    for texto in res["mensajes"]:
        print(texto[0].upper() + texto[1:] + ".")
    if not res["unica"]:
        print("No se llama a numpy.linalg.solve: el sistema no tiene solución única.")
        return
    x = res["solucion"]
    decimales = "" if all(e.is_integer for e in x) else " ≈ " + _aprox(x, n_cifras)
    print("Solución única: " + ", ".join(f"{v} = {e}" for v, e in zip(nombres, x)) + decimales)
    print("Comprobación en las ecuaciones originales (residuo = b - A x):")
    for i in range(m):
        print(f"  ecuación {i + 1}: {_sustitucion(A.row(i), x)} = {(A.row(i) * x)[0]};  "
              f"b{_sub(i + 1)} = {b[i]};  residuo {res['residuo'][i]}")
    print("Residuo nulo: la solución cumple todas las ecuaciones." if res["residuo"] == sp.zeros(m, 1)
          else "El residuo NO es cero: revisa las operaciones.")
    xs = res["numpy"]
    if m != n:
        print(f"numpy.linalg.solve solo acepta sistemas cuadrados: se le dan las ecuaciones "
              f"{_lista(res['ecuaciones_numpy'])}, que son independientes; las demás son combinación de ellas.")
    if xs is None:
        print("Aviso: " + res["aviso_numpy"] + ".")
        return
    diferencia = float(np.max(np.abs(xs - _flotante(x).ravel())))
    print("numpy.linalg.solve (punto flotante): (" + ", ".join(cifras(v, n_cifras) for v in xs) + ")", end="; ")
    if diferencia == 0:
        print("coincide con la solución exacta.")
    elif diferencia < 1e-9 * max(1.0, float(np.max(np.abs(xs)))):
        print(f"coincide con la exacta salvo {cifras(diferencia, 2)}, error de redondeo del punto flotante.")
    else:
        print(f"difiere de la exacta en {cifras(diferencia, 2)}: el sistema está mal condicionado.")


def inversa_gj(A, pasos=False, n_cifras=4):
    """A⁻¹ exacta por Gauss-Jordan sobre [A | I]; ValueError si A no es cuadrada o no tiene inversa."""
    A = matriz_limitada(A, "A")
    if A.rows != A.cols:
        raise ValueError(f"la inversa solo existe para matrices cuadradas (esta es {dims(A)})")
    n = A.rows
    M = A.row_join(sp.eye(n))
    el = _elimina(M, n, con_b=False)
    if pasos:
        muestra_aumentada("[A | I]", M, n)
        _imprime_pasos(el["pasos"], corte=n)
    if el["r"] < n:
        cols = el["libres"]
        cuales = f"la columna {cols[0] + 1}" if len(cols) == 1 else f"las columnas {_lista([c + 1 for c in cols])}"
        raise ValueError(f"A no tiene inversa: al aplicar Gauss-Jordan a [A | I], {cuales} de A se "
                         f"{'queda' if len(cols) == 1 else 'quedan'} sin pivote (rango {el['r']} < {n}), así que la parte "
                         "izquierda nunca llega a I (equivale a det A = 0)")
    Ainv = el["R"][:, n:]
    if pasos:
        muestra_aumentada("[I | A⁻¹]", el["R"], n)
        muestra("A⁻¹", Ainv, n_cifras)
        P = A * Ainv
        muestra("A·A⁻¹", P)
        print("Verificación: A·A⁻¹ = I" if P == sp.eye(n) else "Verificación FALLÓ: A·A⁻¹ ≠ I")
    return Ainv

# %% [markdown]
# **Recorrido paso a paso.** El sistema $x+y+z=6$, $2x+2y+3z=15$, $x+3y+2z=13$. Al limpiar la columna 1, el pivote de la columna 2 queda en cero y hay que intercambiar renglones. Lee cada operación y la matriz que deja; al final, la comprobación en las ecuaciones originales y la comparación con `numpy`.

# %%
_recorrido = gauss_jordan("1, 1, 1, 6; 2, 2, 3, 15; 1, 3, 2, 13")

# %% [markdown]
# La misma eliminación da la inversa: Gauss-Jordan sobre $[A\,|\,I]$ con $A=\begin{bmatrix}1&2\\3&4\end{bmatrix}$. La inversa sale con fracciones exactas.

# %%
inversa_gj("1, 2; 3, 4", pasos=True);

# %% [markdown]
# **Calculadora.** Escribe la matriz aumentada por renglones (la última columna es $\mathbf b$). Prueba también los sistemas de los casos de prueba de abajo, por ejemplo uno con pivote cero (`0, 1, 2; 1, 1, 3`) o uno sin solución (`1, 1, 2; 1, 1, 3`).

# %%
def calculadora_gauss(Ab_txt, ver_pasos, n_cifras):
    with bloque():                 # el resultado aparece completo, de una sola vez
        try:
            Ab = aumentada(Ab_txt)
            res = gauss_jordan(Ab, pasos=False)
        except ValueError as err:
            print("Revisa la entrada:", err); return
        _informe_gauss(Ab, res, n_cifras, con_operaciones=ver_pasos)


widgets.interact(calculadora_gauss,
    Ab_txt=widgets.Textarea(value="2, 1, -1, 8\n-3, -1, 2, -11\n-2, 1, 2, -3", description="[A | b] =",
                            continuous_update=False, layout=widgets.Layout(width="340px", height="120px")),
    ver_pasos=widgets.Checkbox(value=True, description="mostrar cada operación"),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False));

# %% [markdown]
# **Otra calculadora de esta sección: inversa por Gauss-Jordan.** Escribe una matriz cuadrada $A$ (hasta $6\times 6$); la calculadora aplica Gauss-Jordan a $[A\,|\,I]$, muestra cada operación y comprueba que $A\,A^{-1}=I$. Si una columna de $A$ se queda sin pivote, te dice que $A$ no tiene inversa.

# %%
def calculadora_inversa_gj(A_txt, n_cifras):
    with bloque():
        try:
            inversa_gj(A_txt, pasos=True, n_cifras=n_cifras)
        except ValueError as err:     # que A no tenga inversa es un resultado, no un error de captura
            print("Resultado:" if str(err).startswith("A no tiene inversa") else "Revisa la entrada:", err)


widgets.interact(calculadora_inversa_gj,
    A_txt=widgets.Textarea(value="1, 2\n3, 4", description="A =", continuous_update=False,
                           layout=widgets.Layout(width="340px", height="90px")),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False));

# %%
# Casos de prueba de la sección 1 (resultado conocido)
def _sin_solve(fn):
    """(resultado de fn(), True si fn() no llamó a numpy.linalg.solve)."""
    antes = _LLAMADAS_SOLVE[0]
    r = fn()
    return r, _LLAMADAS_SOLVE[0] == antes


_GJ = lambda texto: gauss_jordan(texto, pasos=False)
PRUEBAS_1 = [
    ("1, 1, 1, 6; 2, 2, 3, 15; 1, 3, 2, 13 → (1, 2, 3), con el intercambio R2 ↔ R3 y residuo 0",
     lambda: (lambda r: r["solucion"] == vector("1, 2, 3") and "R2 ↔ R3" in r["operaciones"]
                        and r["residuo"] == sp.zeros(3, 1))(_GJ("1, 1, 1, 6; 2, 2, 3, 15; 1, 3, 2, 13"))),
    ("2, 3, 8; 1, -1, -1 → (1, 2), igual que numpy.linalg.solve",
     lambda: (lambda r: r["solucion"] == vector("1, 2") and np.allclose(r["numpy"], [1, 2]))(_GJ("2, 3, 8; 1, -1, -1"))),
    ("2, 1, -1, 8; -3, -1, 2, -11; -2, 1, 2, -3 → (2, 3, -1)",
     lambda: _GJ("2, 1, -1, 8; -3, -1, 2, -11; -2, 1, 2, -3")["solucion"] == vector("2, 3, -1")),
    ("0.5, 0.25, 1; 0.1, 0.2, 0.5 → (1, 2) exacto: todas las cuentas en fracciones",
     lambda: (lambda r: r["solucion"] == vector("1, 2") and all(e.is_Rational for e in r["escalonada"])
                        and r["residuo"] == sp.zeros(2, 1))(_GJ("0.5, 0.25, 1; 0.1, 0.2, 0.5"))),
    ("0, 1, 2; 1, 1, 3 → intercambio en la columna 1 (R1 ↔ R2) y (1, 2)",
     lambda: (lambda r: r["operaciones"][0] == "R1 ↔ R2" and r["solucion"] == vector("1, 2")
                        and any("el pivote de la columna 1 es cero" in t for _, t, _m in r["pasos"]))(_GJ("0, 1, 2; 1, 1, 3"))),
    ("1, 2, 3; 2, 4, 6 → sin solución única (columna 2 sin pivote) y sin llamar a numpy.linalg.solve",
     lambda: (lambda r, limpio: limpio and not r["unica"] and r["sin_pivote"] == [2] and r["numpy"] is None
                                and "sección 4" in r["mensajes"][0])(*_sin_solve(lambda: _GJ("1, 2, 3; 2, 4, 6")))),
    ("1, 1, 2; 1, 1, 3 → aparece el renglón 0 = 1: no hay solución",
     lambda: (lambda r: r["imposible"] == (2, 1) and r["solucion"] is None
                        and any("0 = 1" in t and "no tiene solución" in t for t in r["mensajes"]))(_GJ("1, 1, 2; 1, 1, 3"))),
    ("inversa_gj: [[2,1],[5,3]] → [[3,-1],[-5,2]]; [[1,2],[2,4]] → no tiene inversa",
     lambda: inversa_gj("2, 1; 5, 3") == matriz("3, -1; -5, 2")
             and _mensaje(lambda: inversa_gj("1, 2; 2, 4")).startswith("A no tiene inversa")),
]
corre_pruebas(PRUEBAS_1, 1)

# %% [markdown]
# **Contrasta (sección 1).** Resuelve a mano $x+2y=7$, $3x-y=7$ con la matriz aumentada: escalona, despeja y comprueba en las dos ecuaciones. Escribe tu solución como texto, primero $x$ y luego $y$.

# %%
mi_solucion = None     # escribe tu solución como texto, primero x y luego y; por ejemplo: "4, -1"

if mi_solucion is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    try:
        mia = vector(str(mi_solucion), "tu solución")
    except ValueError as err:
        print("Revisa tu solución:", err)
    else:
        if mia.rows != 2:          # respuesta incompleta: la referencia todavía no se muestra
            print(f"Tu solución debe tener dos números, primero x y luego y (escribiste {mia.rows}): "
                  "complétala y vuelve a ejecutar la celda.")
        else:
            ref = gauss_jordan("1, 2, 7; 3, -1, 7", pasos=False)["solucion"]
            print(f"La calculadora da: x = {ref[0]}, y = {ref[1]}")
            print(compara_vector(mia, ref, ["x", "y"], "Con R2 ← R2 - 3·R1 el segundo renglón queda solo con y; "
                                                        "despeja y, sustituye en la primera ecuación y comprueba en las dos."))
            print("\nEl procedimiento completo:")
            gauss_jordan("1, 2, 7; 3, -1, 7")
