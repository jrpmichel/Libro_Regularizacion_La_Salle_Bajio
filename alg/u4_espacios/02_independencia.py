# ID: ALG-U4-NB02
# Notebook: alg/u4_espacios.ipynb · sección 4.2 independencia lineal
# Repositorio: alg/u4_espacios/02_independencia.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Independencia lineal
#
# **Definición.** Los vectores $\mathbf v_1,\dots,\mathbf v_k$ son *linealmente independientes* si la única manera de obtener el vector cero,
#
# $$c_1\mathbf v_1+c_2\mathbf v_2+\dots+c_k\mathbf v_k=\mathbf 0,$$
#
# es con **todos** los $c_i=0$. Si hay una combinación con algún $c_i\neq 0$ que da $\mathbf 0$ (una *relación de dependencia*), son *linealmente dependientes*.
#
# - **Dependencia ⇔ alguno sobra.** Son dependientes si y solo si alguno de ellos es combinación lineal de los otros: si $c_j\neq 0$ en la relación, despejas $\mathbf v_j$. Quitar ese vector no cambia el generado.
# - **Prueba con el rango.** Con $A=[\,\mathbf v_1\ \cdots\ \mathbf v_k\,]$, la ecuación de arriba es el sistema homogéneo $A\mathbf c=\mathbf 0$. Por eso
#
#   $$\text{independientes}\iff\operatorname{rango}A=k\iff\operatorname{Nul}(A)=\{\mathbf 0\}.$$
#
#   Cada columna sin pivote es un vector que se escribe con los de las columnas pivote, y cada vector de una base de $\operatorname{Nul}(A)$ es una relación de dependencia.
# - **$n$ vectores en $\mathbb R^n$.** $A$ es cuadrada: son independientes si y solo si $\det A\neq 0$.
# - **Más de $n$ vectores en $\mathbb R^n$** siempre son dependientes: el rango no pasa de $n$, que es menor que $k$.
# - **Un conjunto que contiene al $\mathbf 0$** es dependiente: $1\cdot\mathbf 0=\mathbf 0$ es una relación con un coeficiente distinto de cero.
# - **Dos vectores** son dependientes si y solo si son paralelos (uno es múltiplo del otro; el $\mathbf 0$ cuenta como múltiplo de cualquiera). **Con tres o más no basta revisar de dos en dos:** $(1,2,0)$, $(0,1,3)$ y $(2,5,3)$ no tienen ningún par paralelo y aun así son dependientes, porque el tercero está en el plano de los dos primeros.
#
# La calculadora escribe cada relación con coeficientes enteros sin factor común y con el primero distinto de cero positivo, por ejemplo `2·v1 + 1·v2 - 1·v3 = 0`.

# %%
def independencia(vectores):
    """¿Son linealmente independientes? Diccionario con rango, k (número de vectores), n (componentes),
    independientes (True o False), relaciones (base de Nul(A) escalada a enteros sin factor común: cada una es una
    relación de dependencia), textos de las relaciones ('2·v1 + 1·v2 - 1·v3 = 0'), pivote (columnas pivote, desde
    1), redundantes ({j: [(i, c), ...]}: el vector j, que no es pivote, como combinación exacta de los pivote),
    pares_paralelos (pares (i, j) paralelos o con un cero), det (si k = n; None si no) y las formas de A."""
    A = lista_vectores(vectores)
    n, k = A.shape
    el = _elimina(A, k, con_b=False)
    piv = el["pivotes"]
    relaciones = [_primitivo(d) for d in _direcciones(el["R"], piv, k)]
    redundantes = {j + 1: [(c + 1, el["R"][f, j]) for f, c in enumerate(piv)] for j in el["libres"]}
    pares = [(i + 1, j + 1) for i in range(k) for j in range(i + 1, k)
             if sp.Matrix.hstack(A[:, i], A[:, j]).rank() < 2]
    return {"rango": el["r"], "k": k, "n": n, "independientes": el["r"] == k, "relaciones": relaciones,
            "textos": [texto_relacion(c) for c in relaciones], "pivote": [c + 1 for c in piv],
            "redundantes": redundantes, "pares_paralelos": pares, "det": A.det() if k == n else None,
            "A": A, "escalonada": el["E"], "reducida": el["R"], "pasos": el["pasos"]}


def _texto_redundante(j, pares, n_cifras=None):
    """'v3 = 2·v1 + 1·v2' a partir de [(1, 2), (2, 1)]; con n_cifras, los coeficientes no enteros también en decimales."""
    if not pares:
        return f"v{j} = 0 (es el vector cero)"
    coef = [c for _, c in pares]
    texto = f"v{j} = " + _texto_combinacion(coef, nombres=[f"v{i}" for i, _ in pares])
    if n_cifras is not None and not all(c.is_integer for c in coef):
        texto += f"   (coeficientes ≈ {_aprox(coef, n_cifras)})"
    return texto


def _informe_independencia(res, n_cifras, con_pasos=False):
    A, n, k, r = res["A"], res["n"], res["k"], res["rango"]
    uno = k == 1
    _nombra_vectores(A)
    muestra(_nombre_A(k), A, n_cifras, tex=_tex_A(k))
    if con_pasos:
        _imprime_pasos(_solo_gauss(res["pasos"]), n_cifras=n_cifras)
    muestra("forma escalonada", res["escalonada"], n_cifras)
    print(f"rango(A) = {r},  k = {k} {'vector' if uno else 'vectores'} en {_R(n)}")
    if res["det"] is not None and not uno:
        d = res["det"]
        igual = _igual("det A", d, n_cifras)
        print(f"Como k = n = {n}, A es cuadrada: {igual}" +
              ((" ≠ 0" if "(" not in igual else ", distinto de 0,") + " → independientes." if d != 0 else " → dependientes."))
    if k > n:
        print(f"Hay más vectores ({k}) que componentes ({n}): el rango no pasa de {n}, así que tienen que ser dependientes.")
    ceros = [j + 1 for j in range(k) if A[:, j] == sp.zeros(n, 1)]
    if ceros:
        print(("El vector es el cero: es dependiente (1·0 = 0)." if uno else
               f"El conjunto contiene al vector cero ({_lista([f'v{j}' for j in ceros])}): es dependiente."))
    if res["independientes"]:
        if uno:
            print("rango = k = 1: el vector es LINEALMENTE INDEPENDIENTE (no es el cero); Nul(A) = {0}.")
        else:
            print(f"rango = k = {k}: los vectores son LINEALMENTE INDEPENDIENTES; Nul(A) = {{0}}: la única combinación "
                  "que da 0 es la que tiene todos los coeficientes en 0.")
        return
    print(f"rango = {r} < k = {k}: " + ("el vector es LINEALMENTE DEPENDIENTE." if uno else
                                        "los vectores son LINEALMENTE DEPENDIENTES."))
    piv = res["pivote"]
    if piv:
        print(f"{'Columna pivote' if len(piv) == 1 else 'Columnas pivote'}: {_lista(piv)}; cada vector que no es pivote "
              f"se escribe con {'el vector pivote' if len(piv) == 1 else 'los vectores pivote'}:")
        for j, pares in res["redundantes"].items():
            print(f"  {_texto_redundante(j, pares, n_cifras)}")
    else:
        print("Ninguna columna pivote: " + ("el vector es el cero." if uno else "todos los vectores son el cero."))
    dim = len(res["relaciones"])
    print(f"{'Relación' if dim == 1 else 'Relaciones'} de dependencia (base de Nul(A), dim {dim} = k - rango = {k} - {r}):")
    for c, texto in zip(res["relaciones"], res["textos"]):
        print(f"  {texto}" + ("   (comprobado: A·c = 0)" if A * c == sp.zeros(n, 1) else "   (¡A·c ≠ 0!)"))
    if res["pares_paralelos"]:
        print(f"{'Par paralelo' if len(res['pares_paralelos']) == 1 else 'Pares paralelos'} (uno es múltiplo del otro, "
              "o uno es el cero): " + ", ".join(f"(v{i}, v{j})" for i, j in res["pares_paralelos"]))
    elif k >= 2:
        print("Pares paralelos: ninguno.")
        print("Ningún par es paralelo y aun así son dependientes: con tres o más vectores no basta revisar de dos en dos.")

# %% [markdown]
# **Recorrido paso a paso.** Los vectores $(1,2,0)$, $(0,1,3)$ y $(2,5,3)$. Ningún par es paralelo, pero al escalonar $A$ la tercera columna se queda sin pivote: el rango es 2 y aparece la relación $2\mathbf v_1+\mathbf v_2-\mathbf v_3=\mathbf 0$, es decir, $\mathbf v_3=2\mathbf v_1+\mathbf v_2$. El tercero está en el plano de los dos primeros.

# %%
_informe_independencia(independencia("1, 2, 0; 0, 1, 3; 2, 5, 3"), 4, con_pasos=True)

# %% [markdown]
# **Calculadora.** Escribe los vectores **uno por renglón**. Prueba un conjunto con el vector cero, dos vectores paralelos, más vectores que componentes, y tres vectores de $\mathbb R^3$ sin pares paralelos que sí sean independientes.

# %%
def calculadora_independencia(vectores_txt, n_cifras):
    with bloque():                 # el resultado aparece completo, de una sola vez
        try:
            res = independencia(vectores_txt)
        except ValueError as err:
            print("Revisa la entrada:", err); return
        _informe_independencia(res, n_cifras)


_CONTROLES_2 = dict(
    vectores_txt=widgets.Textarea(value="1, 2, 0\n0, 1, 3\n2, 5, 3", description="vectores:", continuous_update=False,
                                  layout=widgets.Layout(width="340px", height="130px")),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_independencia, **_CONTROLES_2);

# %%
# Casos de prueba de la sección 2 (resultado conocido)
PRUEBAS_2 = [
    ("(2, -1), (-6, 3) → dependientes; relación 3·v1 + 1·v2 = 0; el par (1, 2) es paralelo",
     lambda: (lambda r: not r["independientes"] and r["textos"] == ["3·v1 + 1·v2 = 0"]
                        and r["pares_paralelos"] == [(1, 2)])(independencia("2,-1; -6,3"))),
    ("(1, 0, 1), (0, 1, 1), (1, 1, 0) → independientes, det = -2",
     lambda: (lambda r: r["independientes"] and r["det"] == -2 and r["relaciones"] == [])(independencia("1,0,1; 0,1,1; 1,1,0"))),
    ("(1, 2, 0), (0, 1, 3), (2, 5, 3) → dependientes, relación (2, 1, -1), ningún par paralelo",
     lambda: (lambda r: not r["independientes"] and r["relaciones"] == [vector("2, 1, -1")]
                        and r["textos"] == ["2·v1 + 1·v2 - 1·v3 = 0"] and r["pares_paralelos"] == []
                        and r["redundantes"] == {3: [(1, 2), (2, 1)]})(independencia("1,2,0; 0,1,3; 2,5,3"))),
    ("e1, e2, e3, (1, 1, 1) → dependientes, rango 3, relación (1, 1, 1, -1)",
     lambda: (lambda r: not r["independientes"] and r["rango"] == 3 and r["relaciones"] == [vector("1, 1, 1, -1")])(
         independencia("1,0,0; 0,1,0; 0,0,1; 1,1,1"))),
    ("(1, 2, 3), (0, 0, 0) → dependientes (contiene al cero)",
     lambda: (lambda r: not r["independientes"] and r["rango"] == 1 and (1, 2) in r["pares_paralelos"])(
         independencia("1,2,3; 0,0,0"))),
    ("un solo vector no nulo, (4, -2, 7) → independiente",
     lambda: (lambda r: r["independientes"] and r["rango"] == 1 and r["k"] == 1)(independencia("4,-2,7"))),
    ("corrientes (10, -5, -5), (5, 5, -10), (-5, 10, -5) → dependientes, relación 1·v1 - 1·v2 + 1·v3 = 0",
     lambda: (lambda r: not r["independientes"] and r["relaciones"] == [vector("1, -1, 1")] and r["det"] == 0
                        and r["textos"] == ["1·v1 - 1·v2 + 1·v3 = 0"])(independencia("10,-5,-5; 5,5,-10; -5,10,-5"))),
    ("(1, 0, 0, 1), (0, 1, 0, 1), (0, 0, 1, 1) → independientes, rango 3",
     lambda: (lambda r: r["independientes"] and r["rango"] == 3 and r["det"] is None)(
         independencia("1,0,0,1; 0,1,0,1; 0,0,1,1"))),
    ("prueba de humo: la calculadora corre con sus valores iniciales, sin error y sin «Revisa la entrada»",
     lambda: _humo(calculadora_independencia, _CONTROLES_2)),
]
corre_pruebas(PRUEBAS_2, 2)

# %% [markdown]
# **Contrasta (sección 2).** ¿Los vectores $(1,0,2)$, $(0,1,1)$ y $(2,1,5)$ son linealmente independientes? Escalona a mano $A=[\,\mathbf v_1\ \mathbf v_2\ \mathbf v_3\,]$ (los vectores como columnas), cuenta los pivotes y escribe `mi_rango = ` (un entero) y `mi_respuesta = 'independientes'` o `'dependientes'`.

# %%
mi_rango = None         # un entero, por ejemplo 3
mi_respuesta = None     # 'independientes' o 'dependientes' (también 'LI', 'L.D.', 'no son independientes', ...)

# Respuestas que se aceptan, ya normalizadas (minúsculas, sin acentos ni puntuación). Se compara la frase completa.
_RESPUESTAS = {
    "independientes": ["independientes", "independiente", "linealmente independientes", "linealmente independiente",
                       "son independientes", "es independiente", "son linealmente independientes",
                       "es linealmente independiente", "li", "l i", "no son dependientes", "no son linealmente dependientes",
                       "no dependientes", "no son ld"],
    "dependientes": ["dependientes", "dependiente", "linealmente dependientes", "linealmente dependiente",
                     "son dependientes", "es dependiente", "son linealmente dependientes", "es linealmente dependiente",
                     "ld", "l d", "no son independientes", "no son linealmente independientes", "no independientes",
                     "no son li"],
}


def _respuesta_independencia(texto):
    t = normaliza_texto(texto)
    for clave, frases in _RESPUESTAS.items():
        if t in frases:
            return clave
    raise ValueError(f"«{texto}» no es una de las respuestas que reconoce esta celda: escribe «independientes» "
                     "o «dependientes»")


if mi_rango is None or mi_respuesta is None:     # la referencia se muestra después de tu cálculo
    print("Falta tu cálculo: escribe tu rango y tu respuesta arriba y vuelve a ejecutar la celda.")
else:
    try:
        mio = numero(mi_rango, "tu rango")
        mia = _respuesta_independencia(mi_respuesta)
    except ValueError as err:
        print(f"Revisa tu respuesta: {str(err).rstrip('.')}.")
    else:
        if not (mio.is_integer and 0 <= mio <= 3):     # respuesta imposible: la referencia todavía no se muestra
            print(f"El rango de tres vectores de R³ es un entero de 0 a 3 (escribiste {mio}): revisa tu número "
                  "y vuelve a ejecutar la celda.")
        elif (mio == 3) != (mia == "independientes"):  # incoherente: tampoco se muestra la referencia
            print(f"Tu rango ({mio}) y tu respuesta ({mia}) no son congruentes: con tres vectores, son independientes "
                  "si y solo si el rango es 3. Revisa cuál de los dos está mal y vuelve a ejecutar la celda.")
        else:
            ref = independencia("1, 0, 2; 0, 1, 1; 2, 1, 5")
            clave = "independientes" if ref["independientes"] else "dependientes"
            print(f"La calculadora da: rango {ref['rango']}, {clave}" +
                  (f" ({ref['textos'][0]})" if ref["textos"] else "") + ".")
            print("Rango:", "coincide." if mio == ref["rango"] else
                  "NO coincide: escalona con los vectores como COLUMNAS y cuenta los renglones que no quedan en ceros.")
            print("Respuesta:", "coincide." if mia == clave else
                  "NO coincide: son independientes solo si el rango es igual al número de vectores (3).")
            print("\nEl procedimiento completo:")
            _informe_independencia(ref, 4, con_pasos=True)
