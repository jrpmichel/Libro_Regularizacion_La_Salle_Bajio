# ID: ALG-U5-NB02
# Notebook: alg/u5_ortogonalidad.ipynb · sección 5.2 proyección ortogonal sobre una recta y sobre un subespacio
# Repositorio: alg/u5_ortogonalidad/02_proyeccion.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Proyección ortogonal sobre una recta y sobre un subespacio
#
# **Sobre una recta.** La proyección de $\mathbf x$ sobre la recta por el origen de dirección $\mathbf u\neq\mathbf 0$ es el múltiplo $c\,\mathbf u$ que deja el resto $\mathbf x-c\,\mathbf u$ perpendicular a $\mathbf u$. De $(\mathbf x-c\,\mathbf u)\cdot\mathbf u=0$ sale $c=\dfrac{\mathbf x\cdot\mathbf u}{\mathbf u\cdot\mathbf u}$, es decir,
#
# $$\operatorname{proy}_{\mathbf u}\mathbf x=\frac{\mathbf x\cdot\mathbf u}{\mathbf u\cdot\mathbf u}\,\mathbf u,\qquad\text{con matriz}\qquad P=\frac{\mathbf u\,\mathbf u^T}{\mathbf u^T\mathbf u}.$$
#
# **Descomposición ortogonal.** Si $W$ es un subespacio de $\mathbb R^n$, cada $\mathbf x$ se escribe **de una sola manera** como
#
# $$\mathbf x=\mathbf p+\mathbf r,\qquad \mathbf p\ \text{en}\ W,\qquad \mathbf r\ \text{en}\ W^\perp .$$
#
# $\mathbf p=\operatorname{proy}_W\mathbf x$ es la *proyección ortogonal* de $\mathbf x$ sobre $W$ y $\mathbf r=\mathbf x-\mathbf p$ es el *residuo*.
#
# - **Con una base ortogonal** $\{\mathbf u_1,\dots,\mathbf u_k\}$ de $W$, $\mathbf p$ es la suma de las proyecciones sobre cada vector:
#
#   $$\mathbf p=\frac{\mathbf x\cdot\mathbf u_1}{\mathbf u_1\cdot\mathbf u_1}\,\mathbf u_1+\dots+\frac{\mathbf x\cdot\mathbf u_k}{\mathbf u_k\cdot\mathbf u_k}\,\mathbf u_k .$$
# - **Con una base cualquiera** (las columnas independientes de $A$), $\mathbf p=A\mathbf c$, y pedir que el residuo sea perpendicular a cada columna, $A^T(\mathbf x-A\mathbf c)=\mathbf 0$, da las **ecuaciones normales**
#
#   $$A^TA\,\mathbf c=A^T\mathbf x,\qquad \mathbf p=A\mathbf c,\qquad P=A\,(A^TA)^{-1}A^T .$$
#
#   $A^TA$ es invertible porque las columnas de $A$ son independientes. Si la base trae un vector que es combinación de los anteriores, se quita: no cambia $W$.
#
# **Mejor aproximación.** $\mathbf p$ es el punto de $W$ más cercano a $\mathbf x$. Para cualquier otro $\mathbf w$ de $W$, $\mathbf x-\mathbf w=\mathbf r+(\mathbf p-\mathbf w)$ con $\mathbf p-\mathbf w$ en $W$, y por Pitágoras
#
# $$\|\mathbf x-\mathbf w\|^2=\|\mathbf r\|^2+\|\mathbf p-\mathbf w\|^2\ \ge\ \|\mathbf r\|^2 .$$
#
# Por eso la **distancia de $\mathbf x$ a $W$** es $\|\mathbf r\|$.
#
# **La matriz de proyección** cumple $P^2=P$ (proyectar dos veces es proyectar una) y $P^T=P$ (es simétrica). Si $\mathbf x$ ya está en $W$, $P\mathbf x=\mathbf x$ y $\mathbf r=\mathbf 0$.

# %%
def _x_y_base(x, base):
    """Lee x y la base de W (lista de vectores: cada uno se vuelve una columna de A) y revisa los tamaños."""
    A = lista_vectores(base, "la base de W")
    x = vector(x, "x")
    if x.rows != A.rows:
        raise ValueError(f"x tiene {_componentes(x.rows)} y los vectores de la base tienen {A.rows}: deben coincidir")
    return x, A


def proyeccion_recta(x, u):
    """Proyección de x sobre la recta por el origen de dirección u: diccionario con p = (x·u/u·u)·u, r = x - p,
    coef = x·u/u·u, distancia = ‖r‖ (exacta), P = u uᵀ/(uᵀu), x y u. Si u es el vector cero, ValueError."""
    u = vector(u, "u")
    x = vector(x, "x")
    if x.rows != u.rows:
        raise ValueError(f"x tiene {_componentes(x.rows)} y u tiene {u.rows}: deben tener el mismo número de componentes")
    if _es_nulo(u):
        raise ValueError("u no puede ser el vector cero: el 0 no define ninguna recta")
    coef = _punto(x, u) / _punto(u, u)
    p = coef * u
    r = x - p
    return {"p": p, "r": r, "coef": coef, "distancia": sp.sqrt(_punto(r, r)), "P": u * u.T / _punto(u, u),
            "x": x, "u": u}


def proyeccion_subespacio(x, base):
    """Proyección de x sobre W = gen{v1, ..., vk}, exacta. Diccionario con p, r, coeficientes (los de p en los
    vectores usados), usados (sus números, desde 1), distancia (exacta) y distancia_decimal, metodo ('base ortogonal:
    suma de proyecciones' o 'ecuaciones normales'), ATA, ATx, r_ortogonal (Aᵀr = 0), P = A(AᵀA)⁻¹Aᵀ, avisos, A (los
    vectores usados como columnas) y x. Si un vector es combinación de los anteriores, no se usa y se avisa."""
    x, A0 = _x_y_base(x, base)
    n, k = A0.shape
    piv = _elimina(A0, k, con_b=False)["pivotes"]
    if not piv:
        raise ValueError("todos los vectores de la base son el cero: generan {0} y no hay sobre qué proyectar")
    avisos = []
    for j in range(k):
        if j not in piv:
            avisos.append(f"el vector {j + 1} es el cero y no se usa" if _es_nulo(A0[:, j]) else
                          f"el vector {j + 1} es combinación de los anteriores y no se usa")
    A = A0.extract(list(range(n)), piv)
    ATA, ATx = A.T * A, A.T * x
    ortogonal = all(ATA[i, j] == 0 for i in range(ATA.rows) for j in range(ATA.cols) if i != j)
    if ortogonal:
        c = sp.Matrix([ATx[i] / ATA[i, i] for i in range(A.cols)])
        metodo = "base ortogonal: suma de proyecciones"
    else:
        c = ATA.LUsolve(ATx)
        metodo = "ecuaciones normales"
    p = A * c
    r = x - p
    if _es_nulo(r):
        avisos.append("x ya está en W: su proyección es x mismo y r = 0")
    dist = sp.sqrt(_punto(r, r))
    return {"p": p, "r": r, "coeficientes": list(c), "usados": [j + 1 for j in piv], "distancia": dist,
            "distancia_decimal": float(sp.N(dist, 30)), "metodo": metodo, "ATA": ATA, "ATx": ATx,
            "r_ortogonal": es_perpendicular(r, A0), "P": A * ATA.inv() * A.T, "avisos": avisos, "A": A, "x": x,
            "ortogonal": ortogonal}


def es_perpendicular(r, vectores):
    """True si r es perpendicular a cada uno de los vectores (r·v = 0 para todos), en forma exacta."""
    A = lista_vectores(vectores)
    r = vector(r, "r")
    if r.rows != A.rows:
        raise ValueError(f"r tiene {_componentes(r.rows)} y los vectores tienen {A.rows}: deben coincidir")
    return all(_punto(r, A[:, j]) == 0 for j in range(A.cols))


def suma_de_proyecciones(x, base):
    """La suma de las proyecciones de x sobre cada vector de la base, proy_v1 x + ... + proy_vk x. Diccionario con
    suma, terminos (cada proyección), residuo (x - suma), p (la proyección de verdad) y es_la_proyeccion (True solo
    si la suma coincide con p, lo que se garantiza cuando la base es ortogonal)."""
    x, A = _x_y_base(x, base)
    for j in range(A.cols):
        if _es_nulo(A[:, j]):
            raise ValueError(f"el vector {j + 1} es el cero: no se puede proyectar sobre él")
    terminos = [_punto(x, A[:, j]) / _punto(A[:, j], A[:, j]) * A[:, j] for j in range(A.cols)]
    suma = sum(terminos, sp.zeros(A.rows, 1))
    p = proyeccion_subespacio(x, A)["p"]
    return {"suma": suma, "terminos": terminos, "residuo": x - suma, "p": p, "es_la_proyeccion": suma == p}


def _factor_comun(v):
    """(c, w) con v = c·w, c > 0 y w de enteros sin factor común: (-7/6, 7/6, 7/3) -> (7/6, (-1, 1, 2))."""
    w = _primitivo(v)
    i = next(i for i, e in enumerate(w) if e != 0)
    c = v[i] / w[i]
    return (c, w) if c > 0 else (-c, -w)


def _informe_recta(res, n_cifras):
    x, u = res["x"], res["u"]
    print(f"x = {_tupla(x)},  u = {_tupla(u)}")
    print(f"x·u = {_valor(_punto(x, u), n_cifras)},  u·u = {_valor(_punto(u, u), n_cifras)}  →  "
          f"c = (x·u)/(u·u) = {_valor(res['coef'], n_cifras)}")
    print(f"p = c·u = {_tupla(res['p'], n_cifras)}")
    print(f"r = x - p = {_tupla(res['r'], n_cifras)};  r·u = {_valor(_punto(res['r'], u), n_cifras)}"
          + ("  (perpendicular a u)" if _punto(res["r"], u) == 0 else "  (¡no es perpendicular!)"))
    print(f"distancia de x a la recta = ‖r‖ = {_valor(res['distancia'], n_cifras)}")
    muestra("P = u uᵀ/(uᵀu)", res["P"], n_cifras, tex=r"P=\dfrac{\mathbf u\,\mathbf u^T}{\mathbf u^T\mathbf u}")


def _informe_proyeccion(res, n_cifras, comparar_suma=True):
    A, x, p, r = res["A"], res["x"], res["p"], res["r"]
    k = A.cols
    nombres = [f"v{j}" for j in res["usados"]]
    print(f"x = {_tupla(x)}")
    print("Base de W que se usa: " + ", ".join(f"{nm} = {_tupla(A[:, i])}" for i, nm in enumerate(nombres)) +
          f"  (dim W = {k}: {_subespacio(k, A.rows)})")
    for aviso in res["avisos"]:
        if "no se usa" in aviso:
            print("Aviso:", aviso + ".")
    if res["ortogonal"]:
        if k > 1:
            print(f"La base es ortogonal ({_productos_nulos(nombres)}): AᵀA es diagonal y p es la suma de las "
                  "proyecciones sobre cada vector.")
        muestra("AᵀA", res["ATA"], n_cifras, tex=r"A^TA")
        muestra("Aᵀx", res["ATx"], n_cifras, tex=r"A^T\mathbf x")
        for i, nm in enumerate(nombres):
            print(f"  c{i + 1} = (x·{nm})/({nm}·{nm}) = {_cociente(res['ATx'][i], res['ATA'][i, i], n_cifras)}")
    else:
        malos = [f"{nombres[i]}·{nombres[j]} = {_valor(res['ATA'][i, j], n_cifras)}" for i in range(k)
                 for j in range(i + 1, k) if res["ATA"][i, j] != 0]
        print(f"La base NO es ortogonal ({', '.join(malos)}): se usan las ecuaciones normales AᵀA c = Aᵀx.")
        muestra("AᵀA", res["ATA"], n_cifras, tex=r"A^TA")
        muestra("Aᵀx", res["ATx"], n_cifras, tex=r"A^T\mathbf x")
        print(f"Solución: c = {_tupla(sp.Matrix(res['coeficientes']), n_cifras)}")
    print(f"p = " + " ".join(_termino(c, nm, i == 0) for i, (c, nm) in enumerate(zip(res["coeficientes"], nombres))) +
          f" = {_tupla(p, n_cifras)}")
    print(f"r = x - p = {_tupla(r, n_cifras)}")
    print("Comprobación: Aᵀr = 0, así que r es perpendicular a cada vector de W." if res["r_ortogonal"] else
          "Comprobación FALLÓ: Aᵀr ≠ 0.")
    print(f"Distancia de x a W = ‖r‖ = {_valor(res['distancia'], n_cifras)}")
    if "x ya está en W: su proyección es x mismo y r = 0" in res["avisos"]:
        print("x ya está en W: su proyección es x mismo y r = 0.")
    P = res["P"]
    muestra("P = A(AᵀA)⁻¹Aᵀ", P, n_cifras, tex=r"P=A\,(A^TA)^{-1}A^T")
    print(f"P² = P: {'sí' if P * P == P else 'no'};  Pᵀ = P: {'sí' if P.T == P else 'no'};  "
          f"P·x = p: {'sí' if P * x == p else 'no'}")
    if comparar_suma and not res["ortogonal"]:
        s = suma_de_proyecciones(x, A)
        print(f"Ojo: sumar las proyecciones sobre {_coma(nombres)} daría {_tupla(s['suma'], n_cifras)}, que NO es p: "
              "con una base que no es ortogonal, la suma no es la proyección.")


# --- dibujo ---------------------------------------------------------------------------------------
def _base_ortonormal(V):
    """Base ortonormal (en punto flotante, solo para dibujar) de las columnas de V que en decimales se distinguen
    como independientes: la parte perpendicular de cada columna debe medir más de 1e-9 veces la columna (tolerancia
    RELATIVA, que no depende de la escala de los datos)."""
    Q = []
    for v in V.T:
        nv = np.linalg.norm(v)
        if nv == 0 or not np.isfinite(nv):
            continue
        w = v - sum((q @ v) * q for q in Q)
        w = w - sum((q @ w) * q for q in Q)                # segunda pasada: menos error de redondeo
        if np.linalg.norm(w) > 1e-9 * nv:
            Q.append(w / np.linalg.norm(w))
    return Q


def _corte_con_caja(normal, L):
    """Vértices, en orden, del polígono en que el plano normal·x = 0 corta a la caja [-L, L]³ (para dibujarlo)."""
    esquinas = np.array([(a, b, c) for a in (-L, L) for b in (-L, L) for c in (-L, L)], dtype=float)
    puntos = []
    for i in range(8):
        for j in range(i + 1, 8):
            p, q = esquinas[i], esquinas[j]
            if np.count_nonzero(p != q) != 1:          # solo las 12 aristas de la caja
                continue
            fp, fq = normal @ p, normal @ q
            if fp == fq or fp * fq > 0:
                continue
            puntos.append(p + fp / (fp - fq) * (q - p))
    P = np.unique(np.round(np.array(puntos) / L, 12), axis=0) * L
    centro = P.mean(axis=0)
    u = P[0] - centro
    u /= np.linalg.norm(u)
    w = np.cross(normal / np.linalg.norm(normal), u)
    angulos = np.arctan2((P - centro) @ w, (P - centro) @ u)
    return P[np.argsort(angulos)]


def _marca_recta(p, d1, d2, s):
    """Los tres puntos de la marca de ángulo recto en p, entre las direcciones unitarias d1 y d2, de lado s."""
    return np.array([p + s * d1, p + s * d1 + s * d2, p + s * d2])


def dibuja_proyeccion(x, base, n_cifras=4):
    """Dibuja x, su proyección p sobre W y el residuo r = x - p, con la marca de ángulo recto en p. En R² W es una
    recta; en R³, una recta o un pedazo de plano (sombreado). Para n > 3 no dibuja y lo dice. n_cifras: las cifras
    de la distancia en el título."""
    res = proyeccion_subespacio(x, base)
    A = res["A"]
    n, dim = A.rows, A.cols
    if n > 3 or n == 1:
        print(f"No se dibuja: los vectores son de {_R(n)}" + (" y solo se pueden dibujar los de R² y R³." if n > 3 else
              "; en R¹ no hay nada que proyectar."))
        return res
    if dim == n:
        print(f"No se dibuja: W es todo {_R(n)}, así que p = x y r = 0.")
        return res
    xf, pf, rf = (_flotante(v).ravel() for v in (res["x"], res["p"], res["r"]))
    V = _flotante(A)
    Q = _base_ortonormal(V)
    if len(Q) < dim:
        print("No se dibuja: los vectores de la base son casi dependientes (en decimales no se distinguen de "
              "vectores paralelos), así que el dibujo no sería confiable. Los resultados exactos de arriba sí valen.")
        return res
    L = 1.25 * max(1.0, float(np.abs(np.vstack([xf, pf] + [V[:, j] for j in range(dim)])).max()))
    rn, pn = np.linalg.norm(rf), np.linalg.norm(pf)
    d1 = -pf / pn if pn > 1e-9 * L else Q[0]                  # de p hacia el origen (dentro del triángulo)
    s = 0.06 * L
    if rn > 1e-9 * L:
        s = min(s, 0.4 * rn) if pn <= 1e-9 * L else min(s, 0.4 * rn, 0.4 * pn)
    marca = _marca_recta(pf, d1, rf / rn, s) if rn > 1e-9 * L else None
    corto = lambda v: _tupla(v) if _exacta_legible(v) else "≈ " + _aprox(v, n_cifras)
    d = res["distancia"]
    distancia = f"‖r‖ = {_valor(d, n_cifras)}" if len(_txt(d)) <= 28 else f"‖r‖ ≈ {cifras(d, n_cifras)}"
    titulo = f"Proyección de x = {corto(res['x'])} sobre W: p = {corto(res['p'])}; distancia {distancia}"
    if n == 2:
        fig, ax = plt.subplots(figsize=(6.4, 6.4))
        t = np.array([-2.0, 2.0]) * L
        ax.plot(t * Q[0][0], t * Q[0][1], "-.", color="0.7", lw=1.2, label="W: recta por el origen (gris, delgada)")
        _flecha(ax, (0, 0), xf, _COLORES[0], "-", 2.4)
        _flecha(ax, (0, 0), pf, _COLORES[1], "--", 2.4)
        _flecha(ax, pf, xf, _COLORES[2], ":", 2.4)
        ax.plot(*xf, "o", color="black", ms=7, zorder=5)
        ax.plot(*pf, "s", color="black", ms=6, zorder=5)
        for punta, texto, d in ((xf, "x", (6, 6)), (pf, "p", (8, -14)), ((xf + pf) / 2, "r = x - p", (8, 0))):
            ax.annotate(texto, xy=punta, xytext=d, textcoords="offset points", fontsize=11,
                        bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85))
        if marca is not None:
            ax.plot(*marca.T, "-", color="black", lw=1.2)
        ax.plot([], [], "-", color=_COLORES[0], lw=2.4, label="x (flecha continua, punta redonda)")
        ax.plot([], [], "--", color=_COLORES[1], lw=2.4, label="p = proyección (discontinua, punta cuadrada)")
        ax.plot([], [], ":", color=_COLORES[2], lw=2.4, label="r = x - p (punteada)")
        ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect("equal")
        ax.axhline(0, color="black", lw=0.6); ax.axvline(0, color="black", lw=0.6)
        ax.grid(alpha=0.4); ax.set_xlabel("componente 1"); ax.set_ylabel("componente 2")
        lugar = "lower right"
    else:
        fig = plt.figure(figsize=(7, 6.6))
        ax = fig.add_subplot(projection="3d")
        if dim == 2:
            ax.add_collection3d(Poly3DCollection([_corte_con_caja(np.cross(Q[0], Q[1]), L)], facecolor="0.6",
                                                 alpha=0.22, edgecolor="0.4", linewidth=1.2))
            ax.plot([], [], "s", color="0.6", alpha=0.5, ms=10, label="W: plano por el origen (sombreado)")
        else:
            tope = L / float(np.max(np.abs(Q[0])))
            P = np.outer([-tope, tope], Q[0])
            ax.plot(*P.T, "-.", color="0.7", lw=1.2, label="W: recta por el origen (gris, delgada)")
        ax.plot(*np.array([np.zeros(3), xf]).T, "-", color=_COLORES[0], lw=2.4, label="x (continua, punta redonda)")
        ax.plot(*np.array([np.zeros(3), pf]).T, "--", color=_COLORES[1], lw=2.4,
                label="p = proyección (discontinua, punta cuadrada)")
        ax.plot(*np.array([pf, xf]).T, ":", color=_COLORES[2], lw=2.5, label="r = x - p (punteada)")
        ax.scatter(*xf, marker="o", s=45, color="black", depthshade=False)
        ax.scatter(*pf, marker="s", s=40, color="black", depthshade=False)
        if marca is not None:
            ax.plot(*marca.T, "-", color="black", lw=1.2)
        ax.text(*xf, "  x", fontsize=11)
        ax.text(*pf, "  p", fontsize=11)
        ax.text(*((xf + pf) / 2), "  r", fontsize=11)
        ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_zlim(-L, L)
        ax.set_xlabel("componente 1"); ax.set_ylabel("componente 2"); ax.set_zlabel("componente 3")
        ax.set_box_aspect((1, 1, 1))                         # misma escala en los tres ejes
        ax.view_init(elev=22, azim=-55)
        lugar = "upper left"
    ax.set_title(textwrap.fill(titulo, 60), fontsize=10)
    ax.legend(loc=lugar, fontsize=8)
    plt.show()
    return res

# %% [markdown]
# **Recorrido paso a paso.** (1) $\mathbf x=(2,3,4)$ sobre la recta de $\mathbf u=(1,1,1)$: $\mathbf x\cdot\mathbf u=9$ y $\mathbf u\cdot\mathbf u=3$, así que $\mathbf p=3\,(1,1,1)=(3,3,3)$, $\mathbf r=(-1,0,1)$ y la distancia es $\sqrt2$. (2) $\mathbf x=(1,2,3)$ sobre el plano con base **ortogonal** $\{(1,1,0),(1,-1,1)\}$: basta sumar dos proyecciones, $\mathbf p=(13/6,\,5/6,\,2/3)$, $\mathbf r=\tfrac76(-1,1,2)$ y la distancia es $7/\sqrt6$. (3) El mismo $\mathbf x$ sobre el plano con base $\{(1,0,1),(1,1,0)\}$, que **no** es ortogonal: las ecuaciones normales tienen $A^TA=\begin{bmatrix}2&1\\1&2\end{bmatrix}$ y $A^T\mathbf x=(4,3)$, de donde $\mathbf c=(5/3,\,2/3)$, $\mathbf p=(7/3,\,2/3,\,5/3)$ y $\mathbf r=\tfrac43(-1,1,1)$. Al final, el dibujo de este último caso.

# %%
print("(1) Sobre una recta")
_informe_recta(proyeccion_recta("2, 3, 4", "1, 1, 1"), 4)
print("\n" + "=" * 90 + "\n\n(2) Sobre un plano con base ortogonal")
_p2 = proyeccion_subespacio("1, 2, 3", "1, 1, 0; 1, -1, 1")
_informe_proyeccion(_p2, 4)
print("r = ({})·{}".format(*(lambda c, w: (c, _tupla(w)))(*_factor_comun(_p2["r"]))))
print("\n" + "=" * 90 + "\n\n(3) Sobre un plano con base no ortogonal")
_p3 = proyeccion_subespacio("1, 2, 3", "1, 0, 1; 1, 1, 0")
_informe_proyeccion(_p3, 4, comparar_suma=False)
print("r = ({})·{}".format(*(lambda c, w: (c, _tupla(w)))(*_factor_comun(_p3["r"]))))
dibuja_proyeccion("1, 2, 3", "1, 0, 1; 1, 1, 0");

# %% [markdown]
# **Advertencia: con una base que no es ortogonal, sumar las proyecciones no da la proyección.** Con la base $\{(1,0,1),(1,1,0)\}$ del caso (3) y $\mathbf x=(1,2,3)$, la suma de las proyecciones sobre cada vector es $(7/2,\,3/2,\,2)$: está en el plano (es combinación de los dos vectores), pero **no** es $\mathbf p$. Su residuo, $(-5/2,\,1/2,\,1)$, no es perpendicular a $(1,0,1)$ y es más largo que $\mathbf r$. La suma solo funciona con una base ortogonal; para obtenerla está el proceso de Gram-Schmidt (sección 3).

# %%
_s = suma_de_proyecciones("1, 2, 3", "1, 0, 1; 1, 1, 0")
print(f"proy sobre v1 = {_tupla(_s['terminos'][0])},  proy sobre v2 = {_tupla(_s['terminos'][1])}")
print(f"suma = {_tupla(_s['suma'])}   ¿es la proyección p = {_tupla(_s['p'])}? {'sí' if _s['es_la_proyeccion'] else 'NO'}")
print(f"residuo de la suma = x - suma = {_tupla(_s['residuo'])};  su producto con v1 = (1, 0, 1) es "
      f"{_txt(_punto(_s['residuo'], vector('1, 0, 1')))}, no 0")
print(f"‖x - suma‖² = {_valor(_punto(_s['residuo'], _s['residuo']))}  >  ‖r‖² = {_valor(_punto(_p3['r'], _p3['r']))}: "
      "la suma queda más lejos de x que la proyección.")

# %% [markdown]
# **Calculadora.** Escribe la base de $W$ **un vector por renglón** (cada renglón se vuelve una columna de $A$) y el vector $\mathbf x$. Con una base ortogonal suma las proyecciones; con una base cualquiera resuelve las ecuaciones normales. En $\mathbb R^2$ y $\mathbb R^3$ dibuja $\mathbf x$, $\mathbf p$ y $\mathbf r$. Prueba agregar un vector dependiente (por ejemplo `2, 0, 2`), un $\mathbf x$ que ya esté en el plano, una recta de $\mathbb R^2$ y un subespacio de $\mathbb R^4$.

# %%
def calculadora_proyeccion(base_txt, x_txt, dibujar, n_cifras):
    with bloque():                 # el texto aparece completo, de una sola vez; el dibujo va después
        try:
            res = proyeccion_subespacio(x_txt, base_txt)
        except ValueError as err:
            print("Revisa la entrada:", err); return
        _aviso_coma(base_txt, x_txt)
        _informe_proyeccion(res, n_cifras)
    if dibujar:
        dibuja_proyeccion(res["x"], res["A"], n_cifras)


_CONTROLES_2 = dict(
    base_txt=widgets.Textarea(value="1, 0, 1\n1, 1, 0", description="base de W:", continuous_update=False,
                              layout=widgets.Layout(width="340px", height="120px"), style={"description_width": "initial"}),
    x_txt=widgets.Text(value="1, 2, 3", description="x =", continuous_update=False),
    dibujar=widgets.Checkbox(value=True, description="dibujar (R² y R³)"),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_proyeccion, **_CONTROLES_2);

# %%
# Casos de prueba de la sección 2 (resultado conocido)
_W3 = "1, 0, 1; 1, 1, 0"                       # la base no ortogonal del recorrido
_x3 = vector("1, 2, 3")
PRUEBAS_2 = [
    ("proyeccion_recta((2, 3, 4), (1, 1, 1)) → p = (3, 3, 3), r = (-1, 0, 1), distancia √2",
     lambda: (lambda r: r["p"] == vector("3,3,3") and r["r"] == vector("-1,0,1") and r["distancia"] == sp.sqrt(2)
                        and r["coef"] == 3)(proyeccion_recta("2,3,4", "1,1,1"))),
    ("(1, 2, 3) sobre gen{(1, 1, 0), (1, -1, 1)} (base ortogonal) → p = (13/6, 5/6, 2/3), distancia 7/√6",
     lambda: (lambda r: r["p"] == vector("13/6, 5/6, 2/3") and r["metodo"] == "base ortogonal: suma de proyecciones"
                        and r["r"] == sp.Rational(7, 6) * vector("-1,1,2") and r["distancia"] == 7 / sp.sqrt(6))(
         proyeccion_subespacio("1,2,3", "1,1,0; 1,-1,1"))),
    ("(1, 2, 3) sobre gen{(1, 0, 1), (1, 1, 0)} → c = (5/3, 2/3), p = (7/3, 2/3, 5/3), ecuaciones normales, Aᵀr = 0",
     lambda: (lambda r: r["coeficientes"] == [sp.Rational(5, 3), sp.Rational(2, 3)] and r["p"] == vector("7/3, 2/3, 5/3")
                        and r["metodo"] == "ecuaciones normales" and r["r_ortogonal"]
                        and r["ATA"] == sp.Matrix([[2, 1], [1, 2]]) and r["ATx"] == vector("4, 3"))(
         proyeccion_subespacio("1,2,3", _W3))),
    ("misma base: P² = P, Pᵀ = P y P·x = p",
     lambda: (lambda r: r["P"] * r["P"] == r["P"] and r["P"].T == r["P"] and r["P"] * _x3 == r["p"])(
         proyeccion_subespacio(_x3, _W3))),
    ("suma de proyecciones con la base no ortogonal → (7/2, 3/2, 2), no es la proyección; con la base ortogonal, sí",
     lambda: (lambda s, t: s["suma"] == vector("7/2, 3/2, 2") and not s["es_la_proyeccion"] and t["es_la_proyeccion"])(
         suma_de_proyecciones("1,2,3", _W3), suma_de_proyecciones("1,2,3", "1,1,0; 1,-1,1"))),
    ("(3, 1, 2) ya está en gen{(1, 0, 1), (1, 1, 0)} → p = x, r = 0 y el aviso",
     lambda: (lambda r: _es_nulo(r["r"]) and r["p"] == vector("3,1,2") and r["distancia"] == 0
                        and any("ya está en W" in a for a in r["avisos"]))(proyeccion_subespacio("3,1,2", _W3))),
    ("base con un vector dependiente, (1, 0, 1), (2, 0, 2), (1, 1, 0) → aviso del vector 2 y la misma p",
     lambda: (lambda r: r["p"] == vector("7/3, 2/3, 5/3") and r["usados"] == [1, 3]
                        and "el vector 2 es combinación de los anteriores y no se usa" in r["avisos"])(
         proyeccion_subespacio("1,2,3", "1,0,1; 2,0,2; 1,1,0"))),
    ("en R⁴: (1, 3, 4, 4) sobre gen{(1, 1, 1, 1), (0, 1, 2, 3)} → p = (3/2, 5/2, 7/2, 9/2)",
     lambda: proyeccion_subespacio("1,3,4,4", "1,1,1,1; 0,1,2,3")["p"] == vector("3/2, 5/2, 7/2, 9/2")),
    ("mejor aproximación: tres vectores w de W quedan a distancia ≥ ‖x - p‖ (y mayor si w ≠ p)",
     lambda: (lambda r: all(_punto(_x3 - w, _x3 - w) > _punto(r["r"], r["r"])
                            for w in (2 * vector("1,0,1") + vector("1,1,0"), vector("1,0,1") + vector("1,1,0"),
                                      sp.Rational(3, 2) * vector("1,0,1") + sp.Rational(1, 2) * vector("1,1,0")))
              and _punto(_x3 - r["p"], _x3 - r["p"]) == _punto(r["r"], r["r"]))(proyeccion_subespacio(_x3, _W3))),
    ("proyeccion_recta((1, 2), (0, 0)) → ValueError (u = 0 no define una recta)",
     lambda: "vector cero" in _mensaje(lambda: proyeccion_recta("1,2", "0,0"))),
    ("prueba de humo: la calculadora corre con sus valores iniciales (con y sin dibujo, y con una recta de R²)",
     lambda: _humo(calculadora_proyeccion, _CONTROLES_2)
             and _humo(calculadora_proyeccion, _CONTROLES_2, dibujar=False)
             and _humo(calculadora_proyeccion, _CONTROLES_2, base_txt="2, 1", x_txt="1, 3")),
    ("base casi dependiente (1, 0, 0), (1, 1e-12, 0) con x = (1, 1, 1): no truena, da el resultado exacto "
     "(p = (1, 1, 0), r = (0, 0, 1)) y dice que no se dibuja",
     lambda: (lambda t, e: e is None and "= (1, 1, 0)" in t and "r = x - p = (0, 0, 1)" in t and "No se dibuja" in t
                           and "casi dependientes" in t)(
         *_captura(calculadora_proyeccion, "1, 0, 0\n1, 1e-12, 0", "1, 1, 1", True, 4))),
    ("base con números muy chicos, (1e-15, 0) con x = (1, 1): dibuja la recta sin error y da p = (1, 0)",
     lambda: (lambda t, e: e is None and "= (1, 0)" in t and "No se dibuja" not in t)(
         *_captura(calculadora_proyeccion, "1e-15, 0", "1, 1", True, 4))),
    ("la P de proyeccion_recta es u uᵀ/(uᵀu): con u = (3, 1), P = (1/10)[[9, 3], [3, 1]]",
     lambda: proyeccion_recta("1,2", "3,1")["P"] == sp.Matrix([[9, 3], [3, 1]]) / 10),
    ("es_perpendicular: el residuo del recorrido sí lo es; (1, 0, 0) no es perpendicular a (1, 1, 0)",
     lambda: es_perpendicular(proyeccion_subespacio(_x3, _W3)["r"], _W3)
             and not es_perpendicular("1,0,0", "1,1,0") and es_perpendicular("1,-1,0", "1,1,0")),
    ("un vector cero en la base se avisa como cero, no como combinación: (1, 0, 1), (0, 0, 0), (1, 1, 0)",
     lambda: (lambda r: "el vector 2 es el cero y no se usa" in r["avisos"]
                        and not any("combinación" in a for a in r["avisos"]) and r["usados"] == [1, 3])(
         proyeccion_subespacio("1,2,3", "1,0,1; 0,0,0; 1,1,0"))),
    ("residuo de suma_de_proyecciones con la base no ortogonal: x - suma = (-5/2, 1/2, 1)",
     lambda: suma_de_proyecciones("1,2,3", _W3)["residuo"] == vector("-5/2, 1/2, 1")),
]
corre_pruebas(PRUEBAS_2, 2)

# %% [markdown]
# **Contrasta (sección 2).** Calcula a mano la proyección de $\mathbf x=(0,0,6)$ sobre la recta de $\mathbf u=(1,2,2)$: primero $\mathbf x\cdot\mathbf u$ y $\mathbf u\cdot\mathbf u$, después el coeficiente. Escribe `mi_proyeccion = (a, b, c)`; puedes usar fracciones como texto (`"1/3"`) o decimales.

# %%
mi_proyeccion = None     # tres números, por ejemplo (1, -2, 0.5) o ("1/3", "2/3", 1)

if mi_proyeccion is None:     # la referencia se muestra después de tu cálculo
    print("Falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda.")
else:
    try:
        mia = vector(mi_proyeccion, "tu proyección")
    except ValueError as err:
        print(f"Revisa tu proyección: {str(err).rstrip('.')}.")
    else:
        x_c, u_c = vector("0, 0, 6"), vector("1, 2, 2")          # los datos del enunciado
        mf, xf, uf = (_flotante(v).ravel() for v in (mia, x_c, u_c)) if mia.rows == 3 else (None, None, None)
        if mia.rows != 3:          # respuestas imposibles: la referencia todavía no se muestra
            print(f"x y u son de R³, así que la proyección tiene tres componentes (escribiste {mia.rows}): revisa tu "
                  "respuesta y vuelve a ejecutar la celda.")
        elif np.linalg.norm(mf) > np.linalg.norm(xf) * (1 + 1e-9):
            print(f"La proyección no puede medir más que x (‖x‖ = {cifras(np.linalg.norm(xf), 4)} y tu vector mide "
                  f"{cifras(np.linalg.norm(mf), 4)}). Revisa el coeficiente: se divide entre u·u = ‖u‖², no entre ‖u‖. "
                  "Vuelve a ejecutar la celda.")
        elif np.linalg.norm(mf) < 0.01 * np.linalg.norm(xf):
            print("Tu vector es casi el cero. La proyección sobre la recta es cero solo si x es perpendicular a u: "
                  "calcula x·u, revisa y vuelve a ejecutar la celda.")
        else:
            ref = proyeccion_recta(x_c, u_c)
            pf = _flotante(ref["p"]).ravel()
            print(f"La calculadora da: p = {_tupla(ref['p'], 4)}.")
            paralelo = np.linalg.norm(mf - (mf @ uf) / (uf @ uf) * uf) <= 0.01 * np.linalg.norm(mf)   # relativa
            factor = (mf @ pf) / (pf @ pf)
            if all(abs(float(a - b)) <= 0.01 for a, b in zip(mia, ref["p"])):
                print("Coinciden (con tolerancia de 0.01 por componente).")
            elif paralelo and abs(factor - 1) <= 0.02:
                print("Coinciden salvo redondeo (tu vector difiere de p en menos de 2 %; escribe las fracciones exactas "
                      "si puedes).")
            elif paralelo:
                print(f"NO coinciden: tu vector sí está sobre la recta de u, pero es {cifras(factor, 4)} veces p: "
                      "el coeficiente debe ser (x·u)/(u·u).")
            else:
                print("NO coinciden: la proyección sobre una recta es un múltiplo de u; la tuya no lo es.")
            print("\nEl procedimiento completo:")
            _informe_recta(ref, 4)
