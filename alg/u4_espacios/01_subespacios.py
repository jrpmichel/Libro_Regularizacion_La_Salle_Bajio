# ID: ALG-U4-NB01
# Notebook: alg/u4_espacios.ipynb · sección 4.1 subespacios y combinaciones lineales
# Repositorio: alg/u4_espacios/01_subespacios.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Subespacios de $\mathbb R^n$ y combinaciones lineales
#
# **Combinación lineal.** Con vectores $\mathbf v_1,\dots,\mathbf v_k$ de $\mathbb R^n$ y números $c_1,\dots,c_k$, el vector $c_1\mathbf v_1+c_2\mathbf v_2+\dots+c_k\mathbf v_k$ es una *combinación lineal* de ellos, con *coeficientes* $c_1,\dots,c_k$.
#
# **Generado.** $\operatorname{gen}\{\mathbf v_1,\dots,\mathbf v_k\}$ es el conjunto de **todas** las combinaciones lineales de esos vectores. Preguntar si $\mathbf b$ está en el generado es preguntar si el sistema
#
# $$c_1\mathbf v_1+\dots+c_k\mathbf v_k=\mathbf b\qquad\Longleftrightarrow\qquad A\mathbf c=\mathbf b,\quad A=[\,\mathbf v_1\ \cdots\ \mathbf v_k\,]$$
#
# tiene solución: los vectores son las **columnas** de $A$ y la respuesta la dan los rangos, $\operatorname{rango}A=\operatorname{rango}[A\,|\,\mathbf b]$ (unidad 3).
#
# **Subespacio.** Un subconjunto $H$ de $\mathbb R^n$ es un *subespacio* si cumple las tres condiciones:
#
# 1. $\mathbf 0$ está en $H$;
# 2. $H$ es cerrado bajo la suma: si $\mathbf u$ y $\mathbf w$ están en $H$, $\mathbf u+\mathbf w$ también;
# 3. $H$ es cerrado bajo el producto por escalar: si $\mathbf u$ está en $H$ y $c$ es un número, $c\,\mathbf u$ también.
#
# La condición 1 es la primera prueba: si un conjunto no contiene al origen, no es subespacio.
#
# - **El generado siempre es un subespacio**: la suma de dos combinaciones de los $\mathbf v_i$ y un múltiplo de una combinación son otra vez combinaciones de los $\mathbf v_i$.
# - $\operatorname{Nul}(A)=\{\mathbf x:\ A\mathbf x=\mathbf 0\}$ es un subespacio: si $A\mathbf x=\mathbf 0$ y $A\mathbf y=\mathbf 0$, entonces $A(\mathbf x+\mathbf y)=\mathbf 0$ y $A(c\,\mathbf x)=\mathbf 0$.
# - $\operatorname{Col}(A)$, el generado de las columnas de $A$, es un subespacio (es un generado).
# - Las soluciones de $A\mathbf x=\mathbf b$ con $\mathbf b\neq\mathbf 0$ **no** forman un subespacio: $\mathbf x=\mathbf 0$ da $A\mathbf 0=\mathbf 0\neq\mathbf b$, así que no contienen al origen.
#
# **Cuáles hay.** En $\mathbb R^2$: $\{\mathbf 0\}$, las rectas por el origen y $\mathbb R^2$ completo. En $\mathbb R^3$: $\{\mathbf 0\}$, las rectas por el origen, los planos por el origen y $\mathbb R^3$ completo. Una recta o un plano que no pasa por el origen no es subespacio.
#
# **Describir un generado con ecuaciones.** Un vector $\mathbf w$ está en $\operatorname{gen}\{\mathbf v_1,\dots,\mathbf v_k\}$ si y solo si cumple $\mathbf a\cdot\mathbf w=0$ para todo $\mathbf a$ con $\mathbf a\cdot\mathbf v_i=0$ para cada $i$, es decir, para todo $\mathbf a$ de $\operatorname{Nul}(A^T)$. Por eso las ecuaciones del generado se obtienen de una base de $\operatorname{Nul}(A^T)$: si las coordenadas de $\mathbf w$ son $(x,y,z)$, un plano por el origen en $\mathbb R^3$ tiene una ecuación $a x+b y+c z=0$, y una recta por el origen, dos. Si $\operatorname{rango}A=n$, el generado es todo $\mathbb R^n$ y no hay ecuaciones.
#
# La calculadora escribe cada ecuación con coeficientes enteros sin factor común y con el primero distinto de cero positivo.

# %%
def describe_generado(vectores):
    """Qué es gen{v1, ..., vk}: diccionario con dim (= rango de A), base (los vectores dados que son columnas pivote,
    numerados desde 1), ecuaciones (vectores a tales que a·x = 0 describe el generado: base de Nul(Aᵀ), escalada a
    enteros sin factor común y con la primera entrada distinta de cero positiva), es_todo (True si dim = n), n, k,
    la matriz A (los vectores como columnas) y su forma escalonada."""
    A = lista_vectores(vectores)
    n, k = A.shape
    el = _elimina(A, k, con_b=False)
    ecuaciones = [_primitivo(a) for a in base_nulo(A.T)]
    return {"dim": el["r"], "base": [c + 1 for c in el["pivotes"]], "ecuaciones": ecuaciones, "es_todo": el["r"] == n,
            "n": n, "k": k, "A": A, "escalonada": el["E"]}


def en_generado(vectores, b, pasos=False, n_cifras=4):
    """¿Está b en gen{v1, ..., vk}? Diccionario con es (True o False), r (rango de A), r_aum (rango de [A | b]),
    coeficientes (lista exacta c1, ..., ck con b = c1·v1 + ... + ck·vk; si hay infinitas combinaciones, en función
    de los parámetros s, t, ...), particular (la combinación con los parámetros en 0), parametros, las formas
    escalonada y reducida de [A | b] y, si no está, el renglón '0 = c' que lo prueba. Con pasos=True imprime el
    procedimiento."""
    A = lista_vectores(vectores)
    b = vector(b, "b")
    n, k = A.shape
    if b.rows != n:
        raise ValueError(f"b tiene {_componentes(b.rows)} y los vectores tienen {n}: deben tener el mismo número de componentes")
    el = _elimina(A.row_join(b), k, con_b=True)
    res = {"es": el["imposible"] is None, "r": el["r"], "r_aum": el["r_aum"], "A": A, "b": b,
           "escalonada": el["E"], "reducida": el["R"], "pivote": [c + 1 for c in el["pivotes"]],
           "coeficientes": None, "particular": None, "direcciones": [], "parametros": [], "renglon": None,
           "pasos": el["pasos"]}
    if el["imposible"] is not None:
        f, c = el["imposible"]
        res["renglon"] = f"R{f} de la forma escalonada de [A | b] dice 0 = {c}"
    else:
        x, xp, dirs, ps = solucion_general(el["R"], el["pivotes"], k)
        res.update(coeficientes=list(x), particular=list(xp), direcciones=dirs, parametros=ps)
    if pasos:
        _informe_en_generado(res, n_cifras, con_pasos=True)
    return res


def _nombra_vectores(A):
    print(("Vector: " if A.cols == 1 else "Vectores: ") + ", ".join(f"v{j + 1} = {_tupla(A[:, j])}" for j in range(A.cols)) +
          f"  ({A.cols} {'vector' if A.cols == 1 else 'vectores'} de {_R(A.rows)})")


def _tex_A(k):
    """LaTeX de A = [v1 v2 ... vk] con los vectores en negritas y separados."""
    return r"A=[\," + r"\ ".join(rf"\mathbf v_{{{j}}}" for j in range(1, k + 1)) + r"\,]"


def _nombre_A(k):
    return f"A = [{' '.join(f'v{j}' for j in range(1, k + 1))}]"


def _informe_generado(desc, n_cifras):
    A, n = desc["A"], desc["n"]
    _nombra_vectores(A)
    muestra(_nombre_A(A.cols), A, n_cifras, tex=_tex_A(A.cols))
    print(f"rango(A) = {desc['dim']} → dim gen = {desc['dim']}: el generado es {_subespacio(desc['dim'], n)}.")
    base, sobran = desc["base"], desc["k"] - len(desc["base"])
    if base:
        print(f"{'Columna pivote' if len(base) == 1 else 'Columnas pivote'}: {_lista(base)} → "
              f"{_lista([f'v{j}' for j in base])} {'basta' if len(base) == 1 else 'bastan'} para generarlo"
              + ("." if not sobran else
                 f"; {'el vector que sobra es combinación' if sobran == 1 else 'los que sobran son combinaciones'} "
                 f"de {'él' if len(base) == 1 else 'ellos'}."))
    else:
        print("Ninguna columna pivote: " + ("el vector es el cero." if desc["k"] == 1 else "todos los vectores son el cero."))
    if desc["es_todo"]:
        print(f"Como rango(A) = n = {n}, no hace falta ninguna ecuación: todo vector de {_R(n)} está en el generado.")
    else:
        nombres = nombres_incognitas(n)
        cuantas = len(desc["ecuaciones"])
        print(f"{'Ecuación' if cuantas == 1 else 'Ecuaciones'} del generado ({cuantas} = n - rango = {n} - {desc['dim']}; "
              f"de una base de Nul(Aᵀ)), en las coordenadas ({_coma(nombres)}):")
        for a in desc["ecuaciones"]:
            ceros = all((a.T * A[:, j])[0] == 0 for j in range(A.cols))
            todos = "la cumple v1" if A.cols == 1 else "la cumplen todos los v_i"
            print(f"  {_ecuacion_homogenea(a, nombres)}" + (f"   ({todos})" if ceros else "   (¡un v_i NO la cumple!)"))


def _informe_en_generado(res, n_cifras, con_pasos=False):
    A, b = res["A"], res["b"]
    n, k = A.shape
    if con_pasos:
        _nombra_vectores(A)
    print(f"¿b = {_tupla(b)} está en gen{{{_coma(f'v{j}' for j in range(1, k + 1))}}}? Se busca c con "
          f"{_texto_combinacion(sp.symbols(' '.join(f'c{j}' for j in range(1, k + 1)), seq=True))} = b, "
          "es decir, A c = b:")
    muestra_aumentada("[A | b]", A.row_join(b), k, tex=r"[A\,|\,\mathbf b]", n_cifras=n_cifras)
    if con_pasos:
        _imprime_pasos(_solo_gauss(res["pasos"]), corte=k, n_cifras=n_cifras)
    muestra_aumentada("forma escalonada de [A | b]", res["escalonada"], k,
                      tex=r"\text{forma escalonada de }[A\,|\,\mathbf b]", n_cifras=n_cifras)
    print(f"rango(A) = {res['r']},  rango[A | b] = {res['r_aum']}")
    if not res["es"]:
        print(f"rango[A | b] = rango(A) + 1: b NO está en el generado; {res['renglon']}, que ninguna combinación cumple.")
        return
    print("rango[A | b] = rango(A): b SÍ está en el generado.")
    muestra_aumentada("forma escalonada reducida", res["reducida"], k, n_cifras=n_cifras)
    c, ps = res["coeficientes"], res["parametros"]
    if ps:
        print(f"Hay infinitas combinaciones ({len(ps)} parámetro{'s' if len(ps) > 1 else ''}, {_coma(ps)}): "
              f"b = {_texto_combinacion(c, ps)}")
        print(f"Por ejemplo, con {'todos los parámetros en 0' if len(ps) > 1 else str(ps[0]) + ' = 0'}: "
              f"b = {_texto_combinacion(res['particular'])}")
    else:
        print(f"Una sola combinación: b = {_texto_combinacion(c)}   (coeficientes {_tupla(sp.Matrix(c), n_cifras)})")
    xp = sp.Matrix(res["particular"])
    print(f"Comprobación: A·{_tupla(xp)} = {_tupla(A * xp)}" + (" = b" if A * xp == b else "  ≠ b: ¡revisa!"))


def _base_ortonormal(V):
    """Base ortonormal (en punto flotante, solo para dibujar) de las columnas independientes de V."""
    Q = []
    for v in V.T:
        w = v - sum((q @ v) * q for q in Q)
        if np.linalg.norm(w) > 1e-9 * max(1.0, np.linalg.norm(v)):
            Q.append(w / np.linalg.norm(w))
    return Q


def _proyeccion_exacta(A, b):
    """Proyección ortogonal exacta de b sobre Col(A) (con las columnas pivote): solo para dibujar el segmento que
    separa a b del generado; la decisión de si b está o no la toman los rangos."""
    el = _elimina(A, A.cols, con_b=False)
    if not el["pivotes"]:
        return sp.zeros(A.rows, 1)
    B = A.extract(list(range(A.rows)), el["pivotes"])
    return B * (B.T * B).inv() * B.T * b


def _dentro_caja(p, d, L):
    """Intervalo de t en el que p + t·d queda dentro de la caja [-L, L]³ (vacío si t0 >= t1)."""
    t0, t1 = -np.inf, np.inf
    for k in range(3):
        if abs(d[k]) > 1e-12:
            a, c = sorted(((-L - p[k]) / d[k], (L - p[k]) / d[k]))
            t0, t1 = max(t0, a), min(t1, c)
        elif abs(p[k]) > L:
            return 0.0, 0.0
    return t0, t1


def _corte_con_caja(normal, L):
    """Vértices, en orden, del polígono en que el plano normal·x = 0 corta a la caja [-L, L]³ (para dibujarlo)."""
    esquinas = np.array([(x, y, z) for x in (-L, L) for y in (-L, L) for z in (-L, L)], dtype=float)
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
    P = np.unique(np.round(np.array(puntos), 12), axis=0)
    centro = P.mean(axis=0)
    u = P[0] - centro
    u /= np.linalg.norm(u)
    w = np.cross(normal / np.linalg.norm(normal), u)
    angulos = np.arctan2((P - centro) @ w, (P - centro) @ u)
    return P[np.argsort(angulos)]


def dibuja_generado(vectores, b=None):
    """Dibuja gen{v1, ..., vk}: en R² la recta o el plano entero; en R³ la recta o el pedazo del plano que cabe en
    la caja, con los vectores dados y, si se da b, el punto b (círculo negro si está en el generado, cruz si no).
    Para n > 3 no dibuja y lo dice."""
    desc = describe_generado(vectores)
    A, n, dim = desc["A"], desc["n"], desc["dim"]
    res = en_generado(A, b) if b is not None else None
    if n > 3 or n == 1:
        print(f"No se dibuja: los vectores son de {_R(n)}" + (" y solo se pueden dibujar los de R² y R³." if n > 3 else
              "; en R¹ el generado es {0} o toda la recta real."))
        return desc
    V = _flotante(A)
    bf = _flotante(b if isinstance(b, sp.MatrixBase) else vector(b, "b")).ravel() if b is not None else None
    puntos = [np.zeros(n)] + [V[:, j] for j in range(A.cols)] + ([bf] if bf is not None else [])
    L = 1.25 * max(1.0, float(np.abs(np.vstack(puntos)).max()))
    Q = _base_ortonormal(_flotante(A.extract(list(range(n)), [j - 1 for j in desc["base"]]))) if dim else []
    titulo = f"gen{{{_coma(f'v{j}' for j in range(1, A.cols + 1))}}}: {_subespacio(dim, n)}"
    if dim == n == 3:
        titulo += " (el generado llena todo el espacio: no hay plano ni recta que dibujar)"
    if n == 2:
        fig, ax = plt.subplots(figsize=(5.8, 5.8))
        if dim == 2:
            ax.add_patch(plt.Rectangle((-L, -L), 2 * L, 2 * L, facecolor="0.85", edgecolor="0.5", hatch="//",
                                       alpha=0.5, label="el generado: todo el plano"))
        elif dim == 1:
            t = np.array([-2.0, 2.0]) * L
            ax.plot(t * Q[0][0], t * Q[0][1], "-.", color="0.4", lw=2.2, label="el generado: una recta por el origen")
        else:
            ax.plot(0, 0, "o", color="0.4", ms=10, label="el generado: solo el origen")
        for j in range(A.cols):
            _flecha(ax, (0, 0), V[:, j], _COLORES[j], _ESTILOS[j], 2, f"v{j + 1}")
        if bf is not None:
            if res["es"]:
                ax.plot(*bf, "o", color="black", ms=10, zorder=5, label=f"b = {_tupla(res['b'])}: está")
            else:
                ax.plot(*bf, "X", color="black", ms=12, zorder=5, label=f"b = {_tupla(res['b'])}: NO está")
        ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect("equal")
        ax.axhline(0, color="black", lw=0.6); ax.axvline(0, color="black", lw=0.6)
        ax.grid(alpha=0.4); ax.set_xlabel("x"); ax.set_ylabel("y")
        lugar = "lower right"
    else:
        fig = plt.figure(figsize=(7, 6.5))
        ax = fig.add_subplot(projection="3d")
        if dim == 2:
            ax.add_collection3d(Poly3DCollection([_corte_con_caja(np.cross(Q[0], Q[1]), L)], facecolor="0.6",
                                                 alpha=0.22, edgecolor="0.4", linewidth=1.2))
            rotulo = "el generado: un plano por el origen"
            for c in np.linspace(-1.6 * L, 1.6 * L, 9):                    # rejilla del plano, recortada a la caja
                for d, o in ((Q[0], Q[1]), (Q[1], Q[0])):
                    t0, t1 = _dentro_caja(c * o, d, L)
                    if t0 < t1:
                        P = np.outer([t0, t1], d) + c * o
                        ax.plot(*P.T, "-", color="0.4", lw=0.8, label=rotulo)
                        rotulo = None                                   # una sola entrada en la leyenda
        elif dim == 1:
            t = np.array([-1.0, 1.0]) * L / float(np.max(np.abs(Q[0])))       # la recta, de pared a pared de la caja
            P = np.outer(t, Q[0])
            ax.plot(*P.T, "-.", color="0.4", lw=2.2, label="el generado: una recta por el origen")
        elif dim == 0:
            ax.scatter(0, 0, 0, s=60, color="0.4", depthshade=False, label="el generado: solo el origen")
        for j in range(A.cols):
            ax.plot([0, V[0, j]], [0, V[1, j]], [0, V[2, j]], linestyle=_ESTILOS[j], color=_COLORES[j], lw=2.2,
                    label=f"v{j + 1} = {_tupla(A[:, j])}")
            ax.scatter(*V[:, j], marker=_MARCAS[j], s=45, color=_COLORES[j], depthshade=False)
        if bf is not None:
            if res["es"]:
                ax.scatter(*bf, s=90, marker="o", color="black", depthshade=False, label=f"b = {_tupla(res['b'])}: está")
            else:
                ax.scatter(*bf, s=110, marker="X", color="black", depthshade=False, label=f"b = {_tupla(res['b'])}: NO está")
                if dim > 0:
                    p = _flotante(_proyeccion_exacta(A, res["b"])).ravel()
                    ax.plot(*np.array([bf, p]).T, ":", color="black", lw=1.8, label="distancia de b al generado")
        ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_zlim(-L, L)
        ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
        ax.set_box_aspect((1, 1, 1))                         # misma escala en los tres ejes
        ax.view_init(elev=22, azim=-55)
        lugar = "upper left"
    if res is not None:
        titulo += " · b " + ("está" if res["es"] else "NO está")
    ax.set_title(textwrap.fill(titulo, 60), fontsize=10)
    ax.legend(loc=lugar, fontsize=8)
    plt.show()
    return desc

# %% [markdown]
# **Recorrido paso a paso.** ¿Los vectores $(5,4,1)$ y $(5,4,3)$ están en $\operatorname{gen}\{(1,0,1),(2,1,1)\}$? Se escalona $[A\,|\,\mathbf b]$ con los dos vectores como columnas de $A$. En el primer caso los rangos coinciden y sale la combinación; en el segundo aparece un renglón $0=c$. Al final, la descripción del generado (un plano por el origen y su ecuación) y su dibujo con los dos puntos.

# %%
_r1 = en_generado("1, 0, 1; 2, 1, 1", "5, 4, 1", pasos=True)
print("\n" + "=" * 90 + "\n")
_r2 = en_generado("1, 0, 1; 2, 1, 1", "5, 4, 3", pasos=True)
print("\n" + "=" * 90 + "\n")
_informe_generado(describe_generado("1, 0, 1; 2, 1, 1"), 4)
dibuja_generado("1, 0, 1; 2, 1, 1", "5, 4, 1")
dibuja_generado("1, 0, 1; 2, 1, 1", "5, 4, 3");

# %% [markdown]
# **Calculadora.** Escribe los vectores **uno por renglón** (cada renglón es un vector y se vuelve una columna de $A$) y, si quieres, el vector $\mathbf b$; déjalo vacío para solo describir el generado. Prueba también con vectores de $\mathbb R^2$, por ejemplo `1, 2` y `2, 4` (una recta), o con más vectores que componentes.

# %%
def calculadora_generado(vectores_txt, b_txt, dibujar, n_cifras):
    with bloque():                 # el texto aparece completo, de una sola vez; el dibujo va después
        try:
            A = lista_vectores(vectores_txt)
            b = vector(b_txt, "b") if str(b_txt).strip() else None
            desc = describe_generado(A)
            res = en_generado(A, b) if b is not None else None
        except ValueError as err:
            print("Revisa la entrada:", err); return
        _informe_generado(desc, n_cifras)
        if res is not None:
            print()
            _informe_en_generado(res, n_cifras)
    if dibujar:
        dibuja_generado(A, b)


_CONTROLES_1 = dict(
    vectores_txt=widgets.Textarea(value="1, 0, 1\n2, 1, 1", description="vectores:", continuous_update=False,
                                  layout=widgets.Layout(width="340px", height="120px")),
    b_txt=widgets.Text(value="5, 4, 1", description="b =", continuous_update=False),
    dibujar=widgets.Checkbox(value=True, description="dibujar (R² y R³)"),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_generado, **_CONTROLES_1);

# %%
# Casos de prueba de la sección 1 (resultado conocido)
PRUEBAS_1 = [
    ("v1 = (1, 0, 1), v2 = (2, 1, 1), b = (5, 4, 1) → está, coeficientes (-3, 4), sin parámetros",
     lambda: (lambda r: r["es"] and r["coeficientes"] == [-3, 4] and r["parametros"] == [])(
         en_generado("1,0,1; 2,1,1", "5,4,1"))),
    ("mismos vectores, b = (5, 4, 3) → no está, r = 2, r_aum = 3",
     lambda: (lambda r: not r["es"] and (r["r"], r["r_aum"]) == (2, 3) and r["coeficientes"] is None)(
         en_generado("1,0,1; 2,1,1", "5,4,3"))),
    ("describe_generado((1, 0, 1), (2, 1, 1)) → dim 2, base [1, 2], ecuación (1, -1, -1): x - y - z = 0",
     lambda: (lambda d: d["dim"] == 2 and d["base"] == [1, 2] and d["ecuaciones"] == [vector("1, -1, -1")]
                        and not d["es_todo"])(describe_generado("1,0,1; 2,1,1"))),
    ("describe_generado((1, 2), (2, 4)) → dim 1, base [1], ecuación (2, -1): 2x - y = 0",
     lambda: (lambda d: d["dim"] == 1 and d["base"] == [1] and d["ecuaciones"] == [vector("2, -1")])(
         describe_generado("1,2; 2,4"))),
    ("describe_generado((2, 1, 0), (-1, 0, 1)) → dim 2, ecuación (1, -2, 1)",
     lambda: (lambda d: d["dim"] == 2 and d["ecuaciones"] == [vector("1, -2, 1")])(describe_generado("2,1,0; -1,0,1"))),
    ("describe_generado((1, 0), (1, 1)) → dim 2, es todo R², sin ecuaciones",
     lambda: (lambda d: d["dim"] == 2 and d["es_todo"] and d["ecuaciones"] == [])(describe_generado("1,0; 1,1"))),
    ("b = 0 siempre está: coeficientes (0, 0)",
     lambda: (lambda r: r["es"] and r["coeficientes"] == [0, 0])(en_generado("1,0,1; 2,1,1", "0,0,0"))),
    ("e1, e2, e3 y (1, 1, 1) con b = (2, 3, 4) → está, con 1 parámetro, y toda la familia da b",
     lambda: (lambda r: r["es"] and len(r["parametros"]) == 1
                        and sp.expand(r["A"] * sp.Matrix(r["coeficientes"])) == vector("2, 3, 4"))(
         en_generado("1,0,0; 0,1,0; 0,0,1; 1,1,1", "2,3,4"))),
    ("vectores de distinto tamaño (1, 2) y (1, 2, 3) → ValueError que lo explica",
     lambda: "mismo número de componentes" in _mensaje(lambda: en_generado("1,2; 1,2,3", "1,2"))
             and "vector 2" in _mensaje(lambda: describe_generado("1,2; 1,2,3"))),
    ("prueba de humo: la calculadora corre con sus valores iniciales, sin error y sin «Revisa la entrada»",
     lambda: _humo(calculadora_generado, _CONTROLES_1)),
]
corre_pruebas(PRUEBAS_1, 1)

# %% [markdown]
# **Contrasta (sección 1).** ¿El vector $(1,2,3)$ está en $\operatorname{gen}\{(1,1,1),(0,1,2)\}$? Plantea $c_1(1,1,1)+c_2(0,1,2)=(1,2,3)$ y resuélvelo a mano. Si está, escribe los coeficientes `mis_coeficientes = (c1, c2)`; si no, escribe `mis_coeficientes = 'no'`.

# %%
mis_coeficientes = None     # si está: dos números, por ejemplo (4, -5) o "1/2, 3"; si no está: 'no'

if mis_coeficientes is None:     # la referencia se muestra después de tu cálculo
    print("Falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda.")
else:
    respuesta = None
    if isinstance(mis_coeficientes, str) and normaliza_texto(mis_coeficientes) in (
            "no", "no esta", "no pertenece", "no esta en el generado"):
        respuesta = "no"
    else:
        try:
            respuesta = vector(mis_coeficientes, "tus coeficientes")
        except ValueError as err:
            print(f"Revisa tu respuesta: {str(err).rstrip('.')} (o escribe «no» si el vector no está en el generado).")
        else:
            if respuesta.rows != 2:          # respuesta imposible: la referencia todavía no se muestra
                print(f"Son dos coeficientes, c1 y c2 (escribiste {respuesta.rows}): revisa tu respuesta y vuelve a "
                      "ejecutar la celda.")
                respuesta = None
    if respuesta is not None:
        ref = en_generado("1, 1, 1; 0, 1, 2", "1, 2, 3")
        if ref["es"]:
            print("La calculadora da: sí está, con coeficientes", _tupla(sp.Matrix(ref["coeficientes"])) + ".")
        else:
            print(f"La calculadora da: no está ({ref['renglon']}).")
        if (respuesta == "no") == (not ref["es"]) and (respuesta == "no" or respuesta == sp.Matrix(ref["coeficientes"])):
            print("Coinciden.")
        elif respuesta == "no" or not ref["es"]:
            print("NO coinciden: escalona [A | b] con (1, 1, 1) y (0, 1, 2) como columnas; si no aparece un renglón "
                  "0 = c con c distinto de 0, hay combinación.")
        else:
            print("NO coinciden: la primera componente da c1 directo; sustituye en la segunda para c2 y comprueba "
                  "en la tercera.")
        print("\nEl procedimiento completo:")
        en_generado("1, 1, 1; 0, 1, 2", "1, 2, 3", pasos=True)
