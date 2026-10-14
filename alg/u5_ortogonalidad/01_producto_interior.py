# ID: ALG-U5-NB01
# Notebook: alg/u5_ortogonalidad.ipynb · sección 5.1 producto interior, longitud y ortogonalidad
# Repositorio: alg/u5_ortogonalidad/01_producto_interior.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Producto interior, longitud y ortogonalidad
#
# **Producto interior.** Para $\mathbf u=(u_1,\dots,u_n)$ y $\mathbf v=(v_1,\dots,v_n)$ de $\mathbb R^n$,
#
# $$\mathbf u\cdot\mathbf v=u_1v_1+u_2v_2+\dots+u_nv_n=\mathbf u^T\mathbf v .$$
#
# El resultado es un número, no un vector. Cumple $\mathbf u\cdot\mathbf v=\mathbf v\cdot\mathbf u$, $(c\,\mathbf u)\cdot\mathbf v=c\,(\mathbf u\cdot\mathbf v)$, $\mathbf u\cdot(\mathbf v+\mathbf w)=\mathbf u\cdot\mathbf v+\mathbf u\cdot\mathbf w$ y $\mathbf v\cdot\mathbf v\ge 0$, con $\mathbf v\cdot\mathbf v=0$ solo si $\mathbf v=\mathbf 0$.
#
# **Longitud y distancia.** La *norma* (longitud) es $\|\mathbf v\|=\sqrt{\mathbf v\cdot\mathbf v}$ y la distancia entre $\mathbf u$ y $\mathbf v$ es $\|\mathbf u-\mathbf v\|$. Si $\mathbf v\neq\mathbf 0$, $\mathbf v/\|\mathbf v\|$ es *unitario* (de longitud 1).
#
# **Ángulo.** La desigualdad de Cauchy-Schwarz, $|\mathbf u\cdot\mathbf v|\le\|\mathbf u\|\,\|\mathbf v\|$, garantiza que
#
# $$\cos\theta=\frac{\mathbf u\cdot\mathbf v}{\|\mathbf u\|\,\|\mathbf v\|}$$
#
# está entre $-1$ y $1$ en cualquier $\mathbb R^n$. Por eso define un ángulo $\theta$ entre $0^\circ$ y $180^\circ$ (entre $0$ y $\pi$ radianes) aunque los vectores no se puedan dibujar. Con el vector cero el ángulo no está definido. Ese cociente es la **similitud coseno**: vale 1 si los vectores apuntan igual, 0 si son perpendiculares y $-1$ si apuntan en sentidos opuestos, sin importar su tamaño (por eso se usa para comparar documentos o señales).
#
# **Ortogonalidad.** $\mathbf u$ y $\mathbf v$ son *ortogonales* ($\mathbf u\perp\mathbf v$) si $\mathbf u\cdot\mathbf v=0$. Como $\|\mathbf u+\mathbf v\|^2=\|\mathbf u\|^2+2\,\mathbf u\cdot\mathbf v+\|\mathbf v\|^2$, vale el **teorema de Pitágoras**:
#
# $$\|\mathbf u+\mathbf v\|^2=\|\mathbf u\|^2+\|\mathbf v\|^2\quad\text{si y solo si}\quad\mathbf u\perp\mathbf v .$$
#
# **Conjuntos ortogonales.** $\{\mathbf u_1,\dots,\mathbf u_k\}$ es *ortogonal* si $\mathbf u_i\cdot\mathbf u_j=0$ para cada par $i\neq j$, y *ortonormal* si además cada $\|\mathbf u_i\|=1$. Un conjunto ortogonal **de vectores no nulos** es linealmente independiente: si $c_1\mathbf u_1+\dots+c_k\mathbf u_k=\mathbf 0$, el producto con $\mathbf u_i$ deja $c_i\,(\mathbf u_i\cdot\mathbf u_i)=0$, y como $\mathbf u_i\cdot\mathbf u_i\neq 0$, $c_i=0$. Con los vectores como columnas de $A$, la entrada $(i,j)$ de la *matriz de Gram* $G=A^TA$ es $\mathbf u_i\cdot\mathbf u_j$: el conjunto es ortogonal si $G$ es diagonal y ortonormal si $G=I$.
#
# **Coordenadas en una base ortogonal.** Si $\{\mathbf u_1,\dots,\mathbf u_k\}$ es una base ortogonal de $W$ y $\mathbf x$ está en $W$, no hace falta resolver un sistema: el producto de $\mathbf x=c_1\mathbf u_1+\dots+c_k\mathbf u_k$ con $\mathbf u_i$ deja un solo término, así que
#
# $$c_i=\frac{\mathbf x\cdot\mathbf u_i}{\mathbf u_i\cdot\mathbf u_i}.$$
#
# **Complemento ortogonal.** $W^\perp$ es el conjunto de los vectores perpendiculares a todos los de $W$. Si los vectores que generan $W$ son las columnas de $A$, un $\mathbf y$ es perpendicular a todos ellos si y solo si $A^T\mathbf y=\mathbf 0$, así que
#
# $$W^\perp=\operatorname{Nul}(A^T),\qquad \dim W+\dim W^\perp=n .$$
#
# En $\mathbb R^3$, el complemento de una recta por el origen es el plano perpendicular a ella, y el de un plano por el origen es su recta normal. Los vectores de la base de $W^\perp$ que da la calculadora salen escalados a enteros sin factor común y con la primera entrada distinta de cero positiva; cualquier múltiplo no nulo de cada uno sirve igual (otra herramienta puede darte, por ejemplo, los mismos vectores con el signo cambiado).

# %%
_UNIDADES = {"grados": "grados", "grado": "grados", "radianes": "radianes", "radian": "radianes", "rad": "radianes"}


def _unidad(unidad):
    """'grados' o 'radianes' (acepta mayúsculas, acentos, el singular y 'rad'). No hay unidad por omisión: sin unidad,
    o con otra cosa, ValueError en español."""
    if unidad is None or (isinstance(unidad, str) and not unidad.strip()):
        raise ValueError("indica la unidad: 'grados' o 'radianes'")
    t = normaliza_texto(unidad) if isinstance(unidad, str) else ""
    if t not in _UNIDADES:
        raise ValueError(f"indica la unidad: 'grados' o 'radianes' (escribiste «{unidad}»)")
    return _UNIDADES[t]


def _par(u, v, nu="u", nv="v"):
    """Lee dos vectores y revisa que tengan el mismo número de componentes."""
    u, v = vector(u, nu), vector(v, nv)
    if u.rows != v.rows:
        raise ValueError(f"{nu} tiene {_componentes(u.rows)} y {nv} tiene {v.rows}: deben tener el mismo número "
                         "de componentes")
    return u, v


def producto_interior(u, v):
    """u·v = u1v1 + ... + unvn, exacto."""
    u, v = _par(u, v)
    return _punto(u, v)


def norma(v):
    """‖v‖ = √(v·v), exacta (con raíz si hace falta)."""
    v = vector(v, "el vector")
    return sp.sqrt(_punto(v, v))


def distancia(u, v):
    """Distancia ‖u - v‖, exacta."""
    u, v = _par(u, v)
    return sp.sqrt(_punto(u - v, u - v))


def _coseno(u, v):
    """cos θ = u·v / (‖u‖ ‖v‖), exacto. ValueError si alguno es el vector cero."""
    if _es_nulo(u) or _es_nulo(v):
        raise ValueError("el ángulo no está definido para el vector cero")
    return _punto(u, v) / sp.sqrt(_punto(u, u) * _punto(v, v))


def similitud_coseno(u, v):
    """Similitud coseno u·v / (‖u‖ ‖v‖), exacta: 1 si apuntan igual, 0 si son perpendiculares, -1 si son opuestos.
    Con el vector cero, ValueError «la similitud no está definida para el vector cero»."""
    u, v = _par(u, v)
    if _es_nulo(u) or _es_nulo(v):
        raise ValueError("la similitud no está definida para el vector cero")
    return _coseno(u, v)


def angulo(u, v, unidad=None, n_cifras=4):
    """Ángulo entre u y v. La unidad es obligatoria: angulo(u, v, "grados") o angulo(u, v, "radianes"); sin unidad,
    ValueError «indica la unidad». Devuelve un diccionario con cos (exacto), angulo (exacto si sympy lo simplifica,
    por ejemplo 60 o pi/3; si no, el valor decimal), exacto (True o False), aprox (decimales redondeados a n_cifras),
    unidad y texto ('60 grados', 'π/3 radianes ≈ 1.047', '≈ 73.22 grados')."""
    un = _unidad(unidad)
    u, v = _par(u, v)
    c = _coseno(u, v)
    t = sp.acos(c)
    valor = t if un == "radianes" else t * 180 / sp.pi
    exacto = not valor.has(sp.acos)
    if exacto and not valor.is_Rational:
        valor = sp.simplify(valor)
    aprox = cifras(valor, n_cifras)
    if exacto:
        texto = f"{_txt(valor)} {un}" + ("" if valor.is_integer else f" ≈ {aprox}")
    else:
        texto = f"≈ {aprox} {un}"
    return {"cos": c, "angulo": valor if exacto else float(sp.N(valor, 30)), "exacto": exacto, "aprox": aprox,
            "unidad": un, "texto": texto}


def informe_ortogonalidad(vectores):
    """¿El conjunto es ortogonal? ¿ortonormal? Diccionario con normas (exactas), gram (G = AᵀA exacta), no_ortogonales
    (lista de (i, j, vi·vj) con producto distinto de 0, índices desde 1), ortogonal, ortonormal, hay_cero y A."""
    A = lista_vectores(vectores)
    k = A.cols
    G = A.T * A
    no_ortogonales = [(i + 1, j + 1, G[i, j]) for i in range(k) for j in range(i + 1, k) if G[i, j] != 0]
    hay_cero = any(G[j, j] == 0 for j in range(k))
    ortogonal = not no_ortogonales
    return {"normas": [sp.sqrt(G[j, j]) for j in range(k)], "gram": G, "no_ortogonales": no_ortogonales,
            "ortogonal": ortogonal, "ortonormal": ortogonal and all(G[j, j] == 1 for j in range(k)),
            "hay_cero": hay_cero, "A": A}


def complemento_ortogonal(vectores):
    """W⊥ = Nul(Aᵀ) para W = gen{v1, ..., vk}: diccionario con base (vectores exactos escalados a enteros sin factor
    común, con la primera entrada distinta de cero positiva), dim_W (= rango de A), dim_complemento, n y A."""
    A = lista_vectores(vectores)
    base = [_primitivo(w) for w in base_nulo(A.T)]
    r = _elimina(A, A.cols, con_b=False)["r"]
    return {"base": base, "dim_W": r, "dim_complemento": len(base), "n": A.rows, "A": A}


def coordenadas_ortogonales(base, x):
    """Coordenadas de x en una base ORTOGONAL u1, ..., uk: c_i = x·u_i / (u_i·u_i), exactas, sin resolver ningún
    sistema. Diccionario con coeficientes (lista exacta), en_W (True si x = c1·u1 + ... + ck·uk), suma (ese vector),
    productos (las parejas (x·u_i, u_i·u_i)) y la base como matriz. Si un vector es el cero, o si la base no es
    ortogonal, ValueError."""
    P = lista_vectores(base, "la base")
    x = vector(x, "x")
    n, k = P.shape
    if x.rows != n:
        raise ValueError(f"x tiene {_componentes(x.rows)} y los vectores de la base tienen {n}: deben coincidir")
    ceros = [j + 1 for j in range(k) if _es_nulo(P[:, j])]
    if ceros:
        raise ValueError(f"el vector {ceros[0]} de la base es el cero: una base no puede contener al vector cero")
    G = P.T * P
    malos = [(i + 1, j + 1, G[i, j]) for i in range(k) for j in range(i + 1, k) if G[i, j] != 0]
    if malos:
        i, j, g = malos[0]
        raise ValueError(f"la base no es ortogonal: usa las ecuaciones normales de la sección 2 "
                         f"(v{i}·v{j} = {_valor(g)}, no 0)")
    productos = [(_punto(x, P[:, j]), G[j, j]) for j in range(k)]
    c = [a / b for a, b in productos]
    suma = sum((cj * P[:, j] for j, cj in enumerate(c)), sp.zeros(n, 1))
    return {"coeficientes": c, "en_W": suma == x, "suma": suma, "productos": productos, "base": P, "x": x}


def _informe_ortogonalidad(inf, comp, unidad, n_cifras):
    A, G, k, n = inf["A"], inf["gram"], inf["A"].cols, inf["A"].rows
    _nombra_vectores(A)
    for j, nv in enumerate(inf["normas"], 1):
        print(_igual(f"‖v{j}‖", nv, n_cifras) + ("   (es el vector cero)" if nv == 0 else ""))
    muestra("G = AᵀA", G, n_cifras, tex=r"G=A^TA")
    nombres = [f"v{j}" for j in range(1, k + 1)]
    if k == 1:
        print("G tiene una sola entrada: ‖v1‖². Con un solo vector no hay pares: no hay ángulos que calcular.")
    else:
        print("En la diagonal de G están " + _lista([f"‖{nm}‖²" for nm in nombres]) + "; fuera de ella, los "
              "productos entre vectores distintos (la entrada 1, 2 es v1·v2).")
    for i in range(k):
        for j in range(i + 1, k):
            nombre = f"Ángulo entre v{i + 1} y v{j + 1}"
            if G[i, i] == 0 or G[j, j] == 0:
                print(f"{nombre}: no está definido (hay un vector cero).")
                continue
            a = angulo(A[:, i], A[:, j], unidad, n_cifras)
            c = a["cos"]
            nota = "  (perpendiculares)" if c == 0 else ("  (paralelos, mismo sentido)" if c == 1 else
                                                        ("  (paralelos, sentidos opuestos)" if c == -1 else ""))
            print(f"{nombre}: cos θ = v{i + 1}·v{j + 1}/(‖v{i + 1}‖·‖v{j + 1}‖) = {_valor(c, n_cifras)} → "
                  f"θ = {a['texto']}{nota}")
    if k == 2 and G[0, 0] != 0 and G[1, 1] != 0:
        s = A[:, 0] + A[:, 1]
        izq, der = _punto(s, s), G[0, 0] + G[1, 1]
        if G[0, 1] == 0:
            print(f"Pitágoras: ‖v1 + v2‖² = {_valor(izq, n_cifras)} y ‖v1‖² + ‖v2‖² = "
                  f"{_valor(G[0, 0], n_cifras)} + {_valor(G[1, 1], n_cifras)} = {_valor(der, n_cifras)}: "
                  "se cumple, porque v1·v2 = 0.")
        else:
            print(f"Pitágoras no se cumple: ‖v1 + v2‖² = {_valor(izq, n_cifras)} y ‖v1‖² + ‖v2‖² = "
                  f"{_valor(der, n_cifras)}; la diferencia es 2·v1·v2 = {_valor(2 * G[0, 1], n_cifras)}.")
    if inf["no_ortogonales"]:
        pares = ", ".join(f"v{i}·v{j} = {_valor(g, n_cifras)}" for i, j, g in inf["no_ortogonales"])
        print(f"El conjunto NO es ortogonal: {pares}.")
    elif k == 1:
        print("Un conjunto de un solo vector es ortogonal (no hay pares que revisar).")
    else:
        print(f"El conjunto ES ORTOGONAL: {_productos_nulos(nombres)} (G es diagonal).")
    if inf["ortonormal"]:
        print("Además es ORTONORMAL: " + ("v1 mide 1 (G = [1])." if k == 1 else "cada vector mide 1 (G = I)."))
    elif inf["ortogonal"] and not inf["hay_cero"]:
        print("No es ortonormal: " + ("v1 no mide 1 (v1/‖v1‖ sí)." if k == 1 else
                                     "no todos miden 1. Divide cada vector entre su norma para normalizarlo."))
    if inf["hay_cero"]:
        print("v1 es el vector cero: {0} es dependiente y no sirve de base." if k == 1 else
              "Contiene al vector cero: un conjunto con el cero es dependiente y no sirve de base.")
    elif inf["ortogonal"]:
        print("v1 no es el cero, así que {v1} es linealmente independiente." if k == 1 else
              "Como no contiene al cero, es linealmente independiente.")
    print(f"Complemento ortogonal: W = gen{{{_coma(nombres)}}} tiene dim W = {comp['dim_W']} "
          f"y W⊥ = Nul(Aᵀ) tiene dim W⊥ = {n} - {comp['dim_W']} = {comp['dim_complemento']}.")
    if comp["base"]:
        a_quien = "v1" if k == 1 else _lista(nombres)
        cada = "es perpendicular a" if len(comp["base"]) == 1 else "cada uno es perpendicular a"
        print("Base de W⊥: " + ", ".join(_tupla(w) for w in comp["base"]) +
              (f"   ({cada} {a_quien})" if all(A.T * w == sp.zeros(k, 1) for w in comp["base"])
               else "   (¡alguno NO es perpendicular!)"))
    else:
        print(f"W⊥ = {{0}}: W es todo {_R(n)}.")


def _informe_coordenadas_ortogonales(res, n_cifras):
    P, x = res["base"], res["x"]
    k = P.cols
    print("Base: " + ", ".join(f"u{j} = {_tupla(P[:, j - 1])}" for j in range(1, k + 1)) +
          (f"   (ortogonal: {_productos_nulos([f'u{j}' for j in range(1, k + 1)])})" if k > 1 else ""))
    for j, (a, b) in enumerate(res["productos"], 1):
        print(f"  c{j} = (x·u{j})/(u{j}·u{j}) = {_cociente(a, b, n_cifras)}")
    coef = sp.Matrix(res["coeficientes"])
    if res["en_W"]:
        print(f"Coordenadas: {_tupla(coef, n_cifras)}.  Comprobación: " +
              " ".join(_termino(c, f"u{j}", j == 1) for j, c in enumerate(res["coeficientes"], 1)) +
              f" = {_tupla(res['suma'])} = x")
    else:
        print(f"x NO está en W: la combinación con esos coeficientes da {_tupla(res['suma'], n_cifras)}, que no es x; "
              "es la proyección de x sobre W (sección 2).")

# %% [markdown]
# **Recorrido paso a paso.** Primero, $\mathbf u=(1,-1,2,0)$ y $\mathbf v=(3,1,-1,2)$ de $\mathbb R^4$: $\mathbf u\cdot\mathbf v=3-1-2+0=0$, así que son ortogonales y se cumple Pitágoras: $\|\mathbf u\|^2=6$, $\|\mathbf v\|^2=15$ y $\|\mathbf u+\mathbf v\|^2=21$. Después, el ángulo entre $(1,1,1,1,1,1,1,1,1)$ y $(1,0,0,0,0,0,0,0,0)$, de $\mathbb R^9$, que no se puede dibujar pero sí calcular: $\cos\theta=\dfrac{1}{3\cdot 1}=\tfrac13$ y $\theta\approx 70.53^\circ$ (aquí el ángulo no es uno de los «notables», y la calculadora lo da en decimales). Al final, las coordenadas de $\mathbf x=(3,1,2)$ en la base ortogonal $\{(1,1,1),(1,-1,0),(1,1,-2)\}$, sin resolver ningún sistema: $(2,1,0)$.

# %%
_u, _v = vector("1, -1, 2, 0"), vector("3, 1, -1, 2")
print(f"u = {_tupla(_u)},  v = {_tupla(_v)}")
print(f"u·v = 1·3 + (-1)·1 + 2·(-1) + 0·2 = {producto_interior(_u, _v)}  → u y v son ortogonales")
print(f"‖u‖² = {producto_interior(_u, _u)},  ‖v‖² = {producto_interior(_v, _v)},  u + v = {_tupla(_u + _v)},  "
      f"‖u + v‖² = {producto_interior(_u + _v, _u + _v)}")
print(f"Pitágoras: {producto_interior(_u, _u)} + {producto_interior(_v, _v)} = "
      f"{producto_interior(_u, _u) + producto_interior(_v, _v)} = ‖u + v‖²")
muestra("‖u‖", norma(_u), 4, tex=r"\|\mathbf u\|")
muestra("‖u - v‖", distancia(_u, _v), 4, tex=r"\|\mathbf u-\mathbf v\|")
print()
_unos, _e1 = "1, 1, 1, 1, 1, 1, 1, 1, 1", "1, 0, 0, 0, 0, 0, 0, 0, 0"
_a = angulo(_unos, _e1, "grados")
print(f"Ángulo entre {_tupla(vector(_unos))} y {_tupla(vector(_e1))} (de R⁹): u·v = 1, ‖u‖ = 3, ‖v‖ = 1, "
      f"cos θ = 1/(3·1) = {_txt(_a['cos'])} → θ {_a['texto']} ({angulo(_unos, _e1, 'radianes')['texto']})")
print()
_informe_coordenadas_ortogonales(coordenadas_ortogonales("1, 1, 1; 1, -1, 0; 1, 1, -2", "3, 1, 2"), 4)

# %% [markdown]
# **Calculadora.** Escribe los vectores **uno por renglón** y elige la unidad de los ángulos. La calculadora da las normas, la matriz de Gram $G=A^TA$, el ángulo de cada par, si el conjunto es ortogonal u ortonormal y una base del complemento ortogonal $W^\perp$. Prueba `5/13, 12/13` y `-12/13, 5/13` (ortonormal), un solo vector de $\mathbb R^3$ (su complemento es un plano) y un conjunto con el vector cero.

# %%
def calculadora_ortogonalidad(vectores_txt, unidad, n_cifras):
    with bloque():                 # el resultado aparece completo, de una sola vez
        try:
            inf = informe_ortogonalidad(vectores_txt)
            comp = complemento_ortogonal(inf["A"])
            _unidad(unidad)
        except ValueError as err:
            print("Revisa la entrada:", err); return
        _aviso_coma(vectores_txt)
        _informe_ortogonalidad(inf, comp, unidad, n_cifras)


_CONTROLES_1 = dict(
    vectores_txt=widgets.Textarea(value="1, -1, 2, 0\n3, 1, -1, 2", description="vectores:", continuous_update=False,
                                  layout=widgets.Layout(width="340px", height="120px")),
    unidad=widgets.ToggleButtons(options=["grados", "radianes"], value="grados", description="unidad:"),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_ortogonalidad, **_CONTROLES_1);

# %%
# Casos de prueba de la sección 1 (resultado conocido)
PRUEBAS_1 = [
    ("producto_interior((1, -1, 2, 0), (3, 1, -1, 2)) = 0 y producto_interior((2, -1, 3), (1, 4, 2)) = 4",
     lambda: producto_interior("1,-1,2,0", "3,1,-1,2") == 0 and producto_interior("2,-1,3", "1,4,2") == 4),
    ("norma((1, -1, 2, 0)) = √6 y distancia((1, -1, 2, 0), (3, 1, -1, 2)) = √21",
     lambda: norma("1,-1,2,0") == sp.sqrt(6) and distancia("1,-1,2,0", "3,1,-1,2") == sp.sqrt(21)),
    ("ángulo entre (1, 1, 0) y (-1, 0, 1): 120 grados exacto y 2π/3 en radianes; con (1, ..., 1) y (1, 0, ..., 0) "
     "de R⁹, cos θ = 1/3 y ≈ 70.53 grados (no exacto)",
     lambda: (lambda g, r, d: g["angulo"] == 120 and g["exacto"] and g["cos"] == sp.Rational(-1, 2)
                              and r["angulo"] == 2 * sp.pi / 3 and not d["exacto"] and d["cos"] == sp.Rational(1, 3)
                              and d["aprox"] == "70.53")(
         angulo("1,1,0", "-1,0,1", "grados"), angulo("1,1,0", "-1,0,1", "radianes"),
         angulo("1,1,1,1,1,1,1,1,1", "1,0,0,0,0,0,0,0,0", "grados"))),
    ("ángulo con el vector cero, (0, 0) y (1, 2) → ValueError «no está definido para el vector cero»",
     lambda: "vector cero" in _mensaje(lambda: angulo("0,0", "1,2", "grados"))),
    ("ángulo sin unidad o con otra unidad → ValueError «indica la unidad: 'grados' o 'radianes'»",
     lambda: _mensaje(lambda: angulo("1,0", "0,1")).startswith("indica la unidad: 'grados' o 'radianes'")
             and _mensaje(lambda: angulo("1,0", "0,1", "sexagesimales")).startswith("indica la unidad")),
    ("(1, 1, 1), (1, -1, 0), (1, 1, -2) → ortogonal, no ortonormal, sin cero; G diagonal (3, 2, 6)",
     lambda: (lambda r: r["ortogonal"] and not r["ortonormal"] and not r["hay_cero"]
                        and r["gram"] == sp.diag(3, 2, 6))(informe_ortogonalidad("1,1,1; 1,-1,0; 1,1,-2"))),
    ("(5/13, 12/13), (-12/13, 5/13) → ortonormal (G = I)",
     lambda: (lambda r: r["ortonormal"] and r["gram"] == sp.eye(2))(informe_ortogonalidad("5/13, 12/13; -12/13, 5/13"))),
    ("(1, 2), (2, 1) → no ortogonal, no_ortogonales = [(1, 2, 4)]",
     lambda: (lambda r: not r["ortogonal"] and not r["ortonormal"] and r["no_ortogonales"] == [(1, 2, 4)])(
         informe_ortogonalidad("1,2; 2,1"))),
    ("coordenadas de (3, 1, 2) en la base ortogonal (1, 1, 1), (1, -1, 0), (1, 1, -2) → (2, 1, 0), x en W",
     lambda: (lambda r: r["coeficientes"] == [2, 1, 0] and r["en_W"])(
         coordenadas_ortogonales("1,1,1; 1,-1,0; 1,1,-2", "3,1,2"))),
    ("coordenadas_ortogonales con la base (1, 0), (1, 1), que no es ortogonal → ValueError «la base no es ortogonal»",
     lambda: _mensaje(lambda: coordenadas_ortogonales("1,0; 1,1", "1,1")).startswith("la base no es ortogonal")),
    ("complemento de gen{(1, 1, 1)} → dim 2, cada vector de la base ⊥ (1, 1, 1) y los dos independientes",
     lambda: (lambda r: r["dim_complemento"] == 2 and r["dim_W"] == 1
                        and all(_punto(w, vector("1,1,1")) == 0 for w in r["base"])
                        and sp.Matrix.hstack(*r["base"]).rank() == 2)(complemento_ortogonal("1,1,1"))),
    ("complemento de gen{(1, 2, 0, -1), (0, 1, 1, 2)} → dim 2, perpendicular a los dos, 2 + 2 = 4",
     lambda: (lambda r: r["dim_complemento"] == 2 and r["dim_W"] + r["dim_complemento"] == r["n"] == 4
                        and all(_punto(w, vector(a)) == 0 for w in r["base"] for a in ("1,2,0,-1", "0,1,1,2"))
                        and sp.Matrix.hstack(*r["base"]).rank() == 2)(complemento_ortogonal("1,2,0,-1; 0,1,1,2"))),
    ("similitud_coseno((2, 0), (4, 0)) = 1, ((1, 0), (0, 3)) = 0, ((1, 0), (-2, 0)) = -1; con el cero → ValueError "
     "«la similitud no está definida para el vector cero»",
     lambda: similitud_coseno("2,0", "4,0") == 1 and similitud_coseno("1,0", "0,3") == 0
             and similitud_coseno("1,0", "-2,0") == -1
             and _mensaje(lambda: similitud_coseno("0,0", "1,2")) == "la similitud no está definida para el vector cero"),
    ("prueba de humo: la calculadora corre con sus valores iniciales (en grados y en radianes), sin error y sin "
     "«Revisa la entrada»",
     lambda: _humo(calculadora_ortogonalidad, _CONTROLES_1)
             and _humo(calculadora_ortogonalidad, _CONTROLES_1, unidad="radianes")),
    ("base de W⊥ escalada a enteros y con la primera entrada positiva: gen{(1, 1, 1)} → (1, -1, 0), (1, 0, -1); "
     "gen{(2, 3, 0), (0, 0, 1)} → (3, -2, 0)",
     lambda: complemento_ortogonal("1,1,1")["base"] == [vector("1,-1,0"), vector("1,0,-1")]
             and complemento_ortogonal("2,3,0; 0,0,1")["base"] == [vector("3,-2,0")]),
    ("vectores dependientes (1, 2, 3), (2, 4, 6): dim W = rango = 1 (no 2) y dim W⊥ = 2",
     lambda: (lambda r: r["dim_W"] == 1 and r["dim_complemento"] == 2)(complemento_ortogonal("1,2,3; 2,4,6"))),
    ("x = (1, 2, 3) fuera de W = gen{(1, 1, 0), (1, -1, 0)}: coeficientes (3/2, -1/2), en_W False, suma (1, 2, 0)",
     lambda: (lambda r: r["coeficientes"] == [sp.Rational(3, 2), sp.Rational(-1, 2)] and r["en_W"] is False
                        and r["suma"] == vector("1,2,0"))(coordenadas_ortogonales("1,1,0; 1,-1,0", "1,2,3"))),
    ("coordenadas_ortogonales con un vector cero en la base → ValueError «el vector 2 de la base es el cero»",
     lambda: _mensaje(lambda: coordenadas_ortogonales("1,1,0; 0,0,0", "1,1,0")).startswith(
         "el vector 2 de la base es el cero")),
    ("normas exactas con raíz y vector cero: (1, 2), (2, -1) → normas (√5, √5); (1, 2, 0), (0, 0, 0) → hay_cero",
     lambda: informe_ortogonalidad("1,2; 2,-1")["normas"] == [sp.sqrt(5), sp.sqrt(5)]
             and (lambda r: r["hay_cero"] and r["normas"] == [sp.sqrt(5), 0] and not r["ortonormal"])(
                 informe_ortogonalidad("1,2,0; 0,0,0"))
             and not informe_ortogonalidad("1,1,1; 1,-1,0; 1,1,-2")["hay_cero"]),
]
corre_pruebas(PRUEBAS_1, 1)

# %% [markdown]
# **Contrasta (sección 1).** Calcula a mano el ángulo entre $(1,1,0,0)$ y $(0,1,1,0)$: primero el producto interior y las dos normas, después $\cos\theta$. Escribe `mi_angulo = ` en **grados**.

# %%
mi_angulo = None     # un número en grados, por ejemplo 45 o 37.5

if mi_angulo is None:     # la referencia se muestra después de tu cálculo
    print("Falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda.")
else:
    try:
        mio = numero(mi_angulo, "tu ángulo")
    except ValueError as err:
        print(f"Revisa tu ángulo: {str(err).rstrip('.')}.")
    else:
        if not 0 <= mio <= 180:          # respuesta imposible: la referencia todavía no se muestra
            print(f"El ángulo entre dos vectores está entre 0 y 180 grados (escribiste {_num(mio, 6)}): revisa tu "
                  "cálculo y vuelve a ejecutar la celda.")
        else:
            ref = angulo("1, 1, 0, 0", "0, 1, 1, 0", "grados")
            ref_rad = angulo("1, 1, 0, 0", "0, 1, 1, 0", "radianes")
            print(f"La calculadora da: cos θ = {_txt(ref['cos'])}, θ = {ref['texto']} ({ref_rad['texto']}).")
            if abs(float(mio) - float(ref["angulo"])) <= 0.5:
                print("Coinciden (con tolerancia de 0.5 grados).")
            elif abs(float(mio) - float(ref_rad["angulo"])) <= 0.01:
                print(f"NO coinciden: {_num(mio, 6)} es el ángulo en RADIANES. La celda pide grados: multiplica por "
                      "180/π, o pon tu calculadora en grados al sacar el arco coseno.")
            elif abs(float(mio) - float(ref["cos"])) <= 0.005:
                print(f"NO coinciden: {_num(mio, 6)} es cos θ, no θ. Falta sacar el arco coseno (en grados).")
            else:
                print("NO coinciden: u·v cuenta solo las componentes que son distintas de cero en los dos vectores; "
                      "cada norma es la raíz de la suma de los cuadrados. Después, θ = arccos(u·v/(‖u‖‖v‖)).")
            print("\nEl procedimiento completo:")
            _informe_ortogonalidad(informe_ortogonalidad("1, 1, 0, 0; 0, 1, 1, 0"),
                                   complemento_ortogonal("1, 1, 0, 0; 0, 1, 1, 0"), "grados", 4)
