# ID: ALG-U6-NB04
# Notebook: alg/u6_aplicaciones.ipynb · sección 6.4 armadura por el método de nodos
# Repositorio: alg/u6_aplicaciones/04_armaduras.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 4. Armadura por el método de nodos
#
# **Armadura ideal.** Barras rectas unidas por nodos articulados (pasadores que dejan girar), con las cargas aplicadas solo en los nodos. Así cada barra trabaja solo a lo largo de su eje: o jala o empuja.
#
# **Convención.** La fuerza $F$ de una barra es **positiva si es tensión**: la barra jala a cada uno de sus nodos hacia ella, es decir, hacia el otro extremo. En el nodo $i$, la barra que lo une con el nodo $j$ hace la fuerza $F\,\mathbf u_{ij}$, con $\mathbf u_{ij}=(P_j-P_i)/\|P_j-P_i\|$. Una fuerza negativa es **compresión** (la barra empuja).
#
# **Método de nodos.** Cada nodo está en equilibrio: $\sum F_x=0$ y $\sum F_y=0$. Con $n$ nodos, $b$ barras y $r$ reacciones de los apoyos (un apoyo articulado da dos, $R_x$ y $R_y$; un rodillo da una) hay $2n$ ecuaciones y $b+r$ incógnitas, $A\mathbf x=\mathbf b$, con una columna de $A$ por barra (los vectores unitarios en sus dos nodos) y una por reacción.
#
# - Si $b+r=2n$ y $A$ es invertible, la armadura es **estáticamente determinada**: el equilibrio da todas las fuerzas.
# - Si $b+r<2n$, es un **mecanismo**: le faltan barras o apoyos y se mueve.
# - Si $b+r>2n$, es **hiperestática**: sobran incógnitas, y las ecuaciones de equilibrio no bastan (hace falta cómo se deforman las barras, que verás en mecánica de materiales).
# - **Contar no basta:** con $b+r=2n$, una geometría inestable (por ejemplo, sin nada que impida el movimiento horizontal) da una matriz singular, de rango menor que $2n$.
#
# **Barras de fuerza cero.** En un nodo sin carga ni apoyo al que llegan solo dos barras no colineales, las dos tienen fuerza cero: ninguna puede equilibrar a la otra.
#
# **Cómo se escriben.** Un nodo por renglón como `nombre, x, y` (en metros); una barra como `nombre, nombre`; un apoyo como `nombre, x` o `nombre, y` (un articulado son dos renglones, uno por dirección); una carga como `nombre, Fx, Fy`. Hasta 20 nodos.

# %%
_MAX_NODOS, _MAX_BARRAS = 20, 40
_NOMBRE_NODO = re.compile(r"[A-Za-z][A-Za-z0-9_]{0,7}")


def _nombre_nodo(t, donde):
    t = str(t).strip()
    if not _NOMBRE_NODO.fullmatch(t):
        raise ValueError(f"{donde}: «{t}» no es un nombre de nodo; usa una letra, seguida si quieres de letras o números "
                         "(A, B, N1)")
    return t


def _filas_nombradas(texto, nombre, cosa, ejemplo, maximo):
    if isinstance(texto, dict):
        filas = [[k] + (list(v) if isinstance(v, (list, tuple, np.ndarray, sp.MatrixBase)) else [v])
                 for k, v in texto.items()]
    else:
        filas = _filas_crudas(texto, nombre, cosa, ejemplo)
    if len(filas) > maximo:
        raise ValueError(f"{nombre}: escribiste {len(filas)} renglones y esta calculadora acepta hasta {maximo}")
    return filas


def _lee_nodos(nodos):
    filas = _filas_nombradas(nodos, "los nodos", "el nodo", "A, 0, 0 (nombre, x, y)", _MAX_NODOS)
    N = {}
    for i, f in enumerate(filas, 1):
        if len(f) != 3:
            raise ValueError(f"en los nodos, el renglón {i} debe tener nombre, x, y, por ejemplo: A, 0, 0"
                             f"{_pista_decimal(nodos)}")
        nm = _nombre_nodo(f[0], f"en los nodos, renglón {i}")
        if nm in N:
            raise ValueError(f"el nodo {nm} está escrito dos veces")
        try:
            N[nm] = sp.Matrix([numero(f[1], "x"), numero(f[2], "y")])
        except ValueError as err:
            raise ValueError(f"en el nodo {nm}: {err}") from None
    for a, b in ((a, b) for k, a in enumerate(N) for b in list(N)[k + 1:]):
        if N[a] == N[b]:
            raise ValueError(f"los nodos {a} y {b} están en el mismo punto")
    return N


def _existe(nm, N, donde):
    nm = _nombre_nodo(nm, donde)
    if nm not in N:
        raise ValueError(f"{donde}: el nodo {nm} no está en la lista de nodos ({_coma(N)})")
    return nm


def _lee_barras(barras, N):
    filas = _filas_nombradas(barras, "las barras", "la barra", "A, B (los nombres de sus dos nodos)", _MAX_BARRAS)
    B = []
    for k, f in enumerate(filas, 1):
        if len(f) != 2:
            raise ValueError(f"en las barras, el renglón {k} debe tener los nombres de dos nodos, por ejemplo: A, B")
        i, j = (_existe(e, N, f"en la barra {k}") for e in f)
        if i == j:
            raise ValueError(f"la barra {k} une al nodo {i} consigo mismo")
        if (i, j) in B or (j, i) in B:
            raise ValueError(f"la barra {i}{j} está escrita dos veces")
        B.append((i, j))
    return B


def _lee_apoyos(apoyos, N):
    filas = _filas_nombradas(apoyos, "los apoyos", "el apoyo", "A, x (el nodo y la dirección de la reacción)", 2 * _MAX_NODOS)
    S = []
    for k, f in enumerate(filas, 1):
        if len(f) != 2:
            raise ValueError(f"en los apoyos, el renglón {k} debe tener un nodo y una dirección, por ejemplo: A, x")
        nm = _existe(f[0], N, f"en el apoyo {k}")
        e = normaliza_texto(f[1])
        if e not in ("x", "y"):
            raise ValueError(f"en el apoyo {k}, la dirección debe ser x o y (escribiste «{f[1]}»)")
        if (nm, e) in S:
            raise ValueError(f"el apoyo {nm} en {e} está escrito dos veces")
        S.append((nm, e))
    return S


def _lee_cargas(cargas, N):
    C = {nm: sp.zeros(2, 1) for nm in N}
    if cargas is None or (isinstance(cargas, str) and not cargas.strip()) or (
            isinstance(cargas, (list, tuple, dict)) and len(cargas) == 0):
        return C
    filas = _filas_nombradas(cargas, "las cargas", "la carga", "C, 0, -10 (nodo, Fx, Fy)", 2 * _MAX_NODOS)
    for k, f in enumerate(filas, 1):
        if len(f) != 3:
            raise ValueError(f"en las cargas, el renglón {k} debe tener nodo, Fx, Fy, por ejemplo: C, 0, -10"
                             f"{_pista_decimal(cargas)}")
        nm = _existe(f[0], N, f"en la carga {k}")
        try:
            C[nm] += sp.Matrix([numero(f[1], "Fx"), numero(f[2], "Fy")])     # dos renglones del mismo nodo se suman
        except ValueError as err:
            raise ValueError(f"en la carga {k} (nodo {nm}): {err}") from None
    return C


def _nombre_barra(i, j):
    return f"{i}{j}" if len(i) == len(j) == 1 else f"{i}-{j}"


_EXPLICA = {
    "mecanismo por conteo": "faltan barras o apoyos: la armadura se puede mover (es un mecanismo) y no hay equilibrio para "
                            "una carga cualquiera",
    "hiperestática por conteo": "sobran incógnitas: las ecuaciones de equilibrio no bastan para encontrar las fuerzas "
                                "(hace falta cómo se deforman las barras, que verás en mecánica de materiales)",
    "inestable: la matriz es singular aunque el conteo cuadra": "la geometría permite un movimiento (por ejemplo, nada "
                                                                 "impide que se deslice o gire), así que contar no basta",
    "hiperestática por conteo e inestable (rango < 2n)": "sobran incógnitas y, aun así, la geometría permite un "
                                                         "movimiento: las barras de más no están donde hacen falta",
}


def armadura(nodos, barras, apoyos, cargas, n_cifras=4):
    """Método de nodos, exacto. nodos: dict nombre -> (x, y) o texto 'A, 0, 0' por renglón; barras: lista de pares de
    nombres; apoyos: lista de (nodo, 'x' o 'y'); cargas: dict nodo -> (Fx, Fy) o texto 'C, 0, -10'. Diccionario con
    A (2n × (b + r): vectores unitarios de las barras en sus nodos y una columna por reacción), b (menos las cargas),
    fuerzas (dict barra -> fuerza exacta, positiva = tensión; las longitudes pueden dar raíces), tipo (dict barra ->
    'T', 'C' o '0'), reacciones (dict 'Ax', 'Ay', ... -> valor), longitudes, conteo (n, b, r), rango, diagnostico
    ('determinada', 'mecanismo por conteo', 'hiperestática por conteo', 'hiperestática por conteo e inestable (rango <
    2n)' o 'inestable: la matriz es singular aunque el conteo cuadra') y comprobacion_global (suma de fuerzas externas y momento respecto al primer nodo, que deben ser
    cero) y avisos (por ejemplo, si la armadura está cerca de ser inestable). Si no es determinada, fuerzas, tipo,
    reacciones y comprobacion_global son None."""
    N = _lee_nodos(nodos)
    B, S, C = _lee_barras(barras, N), _lee_apoyos(apoyos, N), _lee_cargas(cargas, N)
    nombres = list(N)
    idx = {nm: k for k, nm in enumerate(nombres)}
    n, b, r = len(N), len(B), len(S)
    M = sp.zeros(2 * n, b + r)                 # columnas de barra multiplicadas por su longitud: entradas racionales
    L = []
    for k, (i, j) in enumerate(B):
        d = N[j] - N[i]
        L.append(sp.sqrt(_punto(d, d)))
        M[2 * idx[i], k], M[2 * idx[i] + 1, k] = d[0], d[1]
        M[2 * idx[j], k], M[2 * idx[j] + 1, k] = -d[0], -d[1]
    for k, (nm, e) in enumerate(S):
        M[2 * idx[nm] + (0 if e == "x" else 1), b + k] = 1
    A = M.copy()
    for k in range(b):
        A[:, k] = M[:, k] / L[k]
    rhs = -sp.Matrix([C[nm][c] for nm in nombres for c in (0, 1)])
    rango = M.rank(iszerofunc=lambda e: e == 0)
    if b + r < 2 * n:
        diag = "mecanismo por conteo"
    elif b + r > 2 * n:
        diag = "hiperestática por conteo" if rango == 2 * n else "hiperestática por conteo e inestable (rango < 2n)"
    elif rango < 2 * n:
        diag = "inestable: la matriz es singular aunque el conteo cuadra"
    else:
        diag = "determinada"
    res = {"A": A, "b": rhs, "fuerzas": None, "tipo": None, "reacciones": None, "comprobacion_global": None,
           "longitudes": {_nombre_barra(i, j): l for (i, j), l in zip(B, L)}, "conteo": (n, b, r), "rango": rango,
           "diagnostico": diag, "nodos": N, "barras": B, "apoyos": S, "cargas": C, "avisos": []}
    if diag != "determinada":
        return res
    cond = _condicion(A)
    if cond > _COND_AVISO:
        res["avisos"] = [f"la armadura está cerca de ser inestable (número de condición de A ≈ {cifras(cond, n_cifras)}): "
                         "las fuerzas pueden salir muy grandes y cambian mucho si la geometría cambia un poco"]
    q = M.LUsolve(rhs)                         # q_k = F_k/L_k en las barras; las reacciones, directas
    fuerzas = {_nombre_barra(i, j): q[k] * L[k] for k, (i, j) in enumerate(B)}
    reacciones = {f"{nm}{e}": q[b + k] for k, (nm, e) in enumerate(S)}
    externas = [(N[nm], C[nm]) for nm in nombres] + [
        (N[nm], sp.Matrix([v, 0]) if e == "x" else sp.Matrix([0, v])) for (nm, e), v in zip(S, reacciones.values())]
    P0 = N[nombres[0]]
    suma = sum((F for _, F in externas), sp.zeros(2, 1))
    mom = sum(((P - P0)[0] * F[1] - (P - P0)[1] * F[0] for P, F in externas), sp.Integer(0))
    res.update(fuerzas=fuerzas, tipo={k: ("T" if v > 0 else "C" if v < 0 else "0") for k, v in fuerzas.items()},
               reacciones=reacciones,
               comprobacion_global={"suma_fuerzas": suma, "momento": mom, "respecto_a": nombres[0],
                                    "cero": _es_nulo(suma) and mom == 0})
    return res


def _informe_armadura(res, unidad, n_cifras):
    N, (n, b, r) = res["nodos"], res["conteo"]
    print("Nodos: " + ", ".join(f"{nm} {_tupla_u(P, '')}" for nm, P in N.items()) + " m")
    print("Barras: " + _coma(_nombre_barra(i, j) for i, j in res["barras"]) +
          ";  apoyos (reacciones): " + (_coma(f"{nm}{e}" for nm, e in res["apoyos"]) or "ninguno"))
    cargadas = [f"{nm} {_tupla_u(F, '')}" for nm, F in res["cargas"].items() if not _es_nulo(F)]
    print(f"Cargas ({unidad}): " + (", ".join(cargadas) if cargadas else "ninguna"))
    print(f"Conteo: n = {n} nodos, b = {b} barras, r = {r} reacciones → 2n = {2 * n} ecuaciones y b + r = {b + r} "
          f"incógnitas;  rango de A = {res['rango']}.")
    diag = res["diagnostico"]
    if diag != "determinada":
        print(f"Diagnóstico: {diag.upper().replace('2N', '2n')} — {_EXPLICA[diag]}.")
        print("No se calculan fuerzas: con esta armadura el equilibrio no da una solución única.")
        return
    print("Diagnóstico: DETERMINADA (b + r = 2n y A es invertible): el equilibrio da todas las fuerzas.")
    if 2 * n <= 8:
        print("Columnas de A: " + _coma([_nombre_barra(i, j) for i, j in res["barras"]] + list(res["reacciones"])) +
              ";  renglones: ΣFx y ΣFy en " + _lista(list(N)) + ".")
        muestra("A", res["A"], n_cifras)
    else:
        print(f"A es de {dims(res['A'])}: no se muestra aquí (está en el resultado de armadura, clave 'A').")
    print(f"Fuerzas en las barras (positiva = tensión, en {unidad}):")
    nombres_tipo = {"T": "tensión", "C": "compresión", "0": "fuerza cero"}
    for k, F in res["fuerzas"].items():
        t = res["tipo"][k]
        print(f"  {k}: {_magnitud('L', res['longitudes'][k], 'm', n_cifras)},  "
              + _magnitud("F", F, unidad, n_cifras) + f"  ({t}, {nombres_tipo[t]})")
    print("Reacciones: " + ",  ".join(_magnitud(k, v, unidad, n_cifras) for k, v in res["reacciones"].items()))
    cg = res["comprobacion_global"]
    print(f"Comprobación global (cargas y reacciones): suma de fuerzas = {_tupla_u(cg['suma_fuerzas'], unidad)}, momento "
          f"respecto a {cg['respecto_a']} = {_ntxt(cg['momento'])} {unidad}·m" + ("  (cero: la armadura completa está en "
                                                                                 "equilibrio)" if cg["cero"] else "  (¡no es cero!)"))
    for aviso in res["avisos"]:
        print("Aviso:", aviso + ".")
    ceros = [k for k, t in res["tipo"].items() if t == "0"]
    if ceros:
        print(f"Barras de fuerza cero: {_coma(ceros)}. No cargan nada con esta carga, pero no se pueden quitar: "
              "rigidizan la armadura y otra carga sí las haría trabajar.")


def dibuja_armadura(resultado, unidad="kN", n_cifras=4):
    """Dibuja la armadura: tensión con línea continua, compresión discontinua y fuerza cero punteada, con su valor y
    T/C; las cargas (flechas negras) y las reacciones (flechas de rayas y puntos). Con más de 12 barras, cada barra
    lleva solo su nombre, y con más de 24, ningún rótulo (los nombres y las fuerzas están en el texto). Si no es
    determinada, dibuja solo la geometría y los apoyos."""
    N, B = resultado["nodos"], resultado["barras"]
    P = {nm: _flotante(v).ravel() for nm, v in N.items()}
    xy = np.array(list(P.values()))
    s = max(float((xy.max(axis=0) - xy.min(axis=0)).max()), 1.0)
    muchas = len(B) > 12                       # con muchas barras, los rótulos con valor se enciman
    sin_rotulos = len(B) > 24                  # y con más de 24, también los nombres
    a = (0.1 if muchas else 0.18) * s
    fig, ax = plt.subplots(figsize=(11.0, 8.6) if muchas else (8.4, 7.0))
    fuerzas, tipo = resultado["fuerzas"], resultado["tipo"]
    estilos = {"T": ("-", _COLORES[0], 2.8), "C": ("--", _COLORES[1], 2.8), "0": (":", "0.45", 2.2)}
    letra = 7.5 if muchas else 8
    for i, j in B:
        k = _nombre_barra(i, j)
        ls, color, lw = estilos[tipo[k]] if tipo else ("-", "0.55", 2.0)
        ax.plot(*np.array([P[i], P[j]]).T, linestyle=ls, color=color, lw=lw, zorder=2)
        if sin_rotulos:
            continue
        texto = k if (not fuerzas or muchas) else f"{k} = {_dec_o_cifras(fuerzas[k], n_cifras)} {unidad} ({tipo[k]})"
        d = P[j] - P[i]
        giro = math.degrees(math.atan2(d[1], d[0]))
        giro = giro - 180 if giro > 90 else (giro + 180 if giro < -90 else giro)      # el texto nunca de cabeza
        ax.annotate(texto, xy=(P[i] + P[j]) / 2, ha="center", va="center", fontsize=letra, zorder=4, rotation=giro,
                    rotation_mode="anchor",
                    bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85))
    for nm, p in P.items():
        ax.plot(*p, "o", color="black", ms=7, zorder=5)
        ax.annotate(nm, xy=p, xytext=(-12, 8), textcoords="offset points", fontsize=8 if len(P) > 10 else 11,
                    fontweight="bold")
    for nm, F in resultado["cargas"].items():
        f = _flotante(F).ravel()
        if np.linalg.norm(f) > 0:
            u = f / np.linalg.norm(f)
            _flecha(ax, P[nm] - a * u, P[nm], "black", "-", 2.6)
            dx, dy = -14 * u                         # el rótulo, más allá de la cola de la flecha (lejos de las barras)
            ax.annotate(f"{_dec_o_cifras(sp.sqrt(_punto(F, F)), n_cifras)} {unidad}", xy=P[nm] - a * u, xytext=(dx, dy),
                        textcoords="offset points", fontsize=9, zorder=6,
                        ha="center" if abs(dx) < 5 else ("right" if dx < 0 else "left"),
                        va="bottom" if dy > 5 else ("top" if dy < -5 else "center"),
                        bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85))
    reac = resultado["reacciones"]
    for (nm, e) in resultado["apoyos"]:
        eje = np.array([1.0, 0.0]) if e == "x" else np.array([0.0, 1.0])
        ax.plot(*(P[nm] - 0.06 * s * eje), ">" if e == "x" else "^", color="0.3", ms=11, zorder=3)
        if reac is not None:
            v = reac[f"{nm}{e}"]
            base = P[nm] - 0.06 * s * eje            # junto al símbolo del apoyo
            cola = base - a * eje
            if v > 0:                                # la flecha apunta en el sentido de la reacción
                _flecha(ax, cola, base, _COLORES[2], "-.", 2.2)
            elif v < 0:
                _flecha(ax, base, cola, _COLORES[2], "-.", 2.2)
            # una reacción cero no se dibuja como flecha: solo su rótulo
            ax.annotate(f"{nm}{e} = {_dec_o_cifras(v, n_cifras)} {unidad}", xy=cola if e == "y" else (cola + P[nm]) / 2,
                        xytext=(0, -10) if e == "y" else (0, 9), textcoords="offset points", fontsize=9,
                        ha="center", va="top" if e == "y" else "bottom",
                        bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85))
    for ls, color, lw, etiqueta in (("-", _COLORES[0], 2.8, "tensión (T): continua"),
                                    ("--", _COLORES[1], 2.8, "compresión (C): discontinua"),
                                    (":", "0.45", 2.2, "fuerza cero: punteada"),
                                    ("-", "black", 2.6, "carga (flecha negra)"),
                                    ("-.", _COLORES[2], 2.2, "reacción (rayas y puntos)")):
        ax.plot([], [], linestyle=ls, color=color, lw=lw, label=etiqueta)
    ax.set_xlim(xy[:, 0].min() - 0.35 * s, xy[:, 0].max() + 0.35 * s)
    ax.set_ylim(xy[:, 1].min() - 0.35 * s, xy[:, 1].max() + 0.35 * s)
    ax.set_aspect("equal"); ax.grid(alpha=0.3)
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)")
    nota = ("" if not muchas else " (los nombres y las fuerzas de las barras están en el texto)" if sin_rotulos else
            " (las fuerzas de cada barra están en el texto)" if fuerzas else "")
    ax.set_title(f"Armadura: {resultado['diagnostico']}" + nota, fontsize=10)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3)
    plt.tight_layout()
    plt.show()

# %% [markdown]
# **Recorrido paso a paso.** (1) Un triángulo: $A(0,0)$ articulado, $B(4,0)$ sobre un rodillo (solo reacción vertical) y $C(2,\,1.5)$ con una carga de 10 kN hacia abajo. Las barras $AC$ y $BC$ miden 2.5 m. En $C$: $\sum F_x=0$ da $F_{AC}=F_{BC}$ y $\sum F_y=0$ da $-0.6\,F_{AC}-0.6\,F_{BC}-10=0$, así que $F_{AC}=F_{BC}=-\tfrac{25}{3}$ kN (compresión); en $B$, $F_{AB}=\tfrac{20}{3}$ kN (tensión), y las reacciones son $A_x=0$, $A_y=B_y=5$ kN. (2) Una armadura Warren: $A(0,0)$, $B(2,0)$, $C(4,0)$, $D(1,2)$, $E(3,2)$, barras $AB$, $BC$, $AD$, $BD$, $BE$, $CE$, $DE$, $A$ articulado, $C$ sobre un rodillo y 6 kN hacia abajo en $B$: $F_{AB}=F_{BC}=1.5$, $F_{AD}=F_{CE}=-1.5\sqrt5$, $F_{BD}=F_{BE}=1.5\sqrt5$, $F_{DE}=-3$ kN, y $A_y=C_y=3$ kN. Las dos cumplen $b+r=2n$ y su matriz es invertible.

# %%
print("(1) Triángulo")
_tri = armadura({"A": (0, 0), "B": (4, 0), "C": (2, 1.5)}, [("A", "B"), ("A", "C"), ("B", "C")],
                [("A", "x"), ("A", "y"), ("B", "y")], {"C": (0, -10)})
_informe_armadura(_tri, "kN", 4)
dibuja_armadura(_tri, "kN")
print("=" * 90 + "\n\n(2) Armadura Warren")
_NODOS_W = "A, 0, 0\nB, 2, 0\nC, 4, 0\nD, 1, 2\nE, 3, 2"
_BARRAS_W = "A, B\nB, C\nA, D\nB, D\nB, E\nC, E\nD, E"
_APOYOS_W = "A, x\nA, y\nC, y"
_warren = armadura(_NODOS_W, _BARRAS_W, _APOYOS_W, "B, 0, -6")
_informe_armadura(_warren, "kN", 4)
dibuja_armadura(_warren, "kN")

# %% [markdown]
# **Calculadora.** Escribe los nodos (`nombre, x, y` en metros), las barras (`nombre, nombre`), los apoyos (`nombre, x` o `nombre, y`; un articulado son dos renglones) y las cargas (`nombre, Fx, Fy`), uno por renglón, y la unidad de fuerza. Hasta 20 nodos. Prueba quitar el apoyo `A, x` (mecanismo), cambiarlo por `C, y` (el conteo cuadra, pero la matriz es singular) y agregar un nodo `D, 4, 1.5` con las barras `C, D` y `B, D`: a $D$ no llega ninguna carga y solo llegan dos barras no colineales, así que las dos quedan con fuerza cero.

# %%
def calculadora_armadura(nodos_txt, barras_txt, apoyos_txt, cargas_txt, unidad, n_cifras):
    with bloque():                 # el texto aparece completo, de una sola vez; el dibujo va después
        try:
            u = nombre_unidad(unidad, "fuerza", "kN o N")
            res = armadura(nodos_txt, barras_txt, apoyos_txt, cargas_txt, n_cifras=n_cifras)
        except ValueError as err:
            _reporta(err); return
        _aviso_coma(nodos_txt, cargas_txt)
        _informe_armadura(res, u, n_cifras)
    dibuja_armadura(res, u, n_cifras)


_caja = lambda alto: widgets.Layout(width="250px", height=alto)
_CONTROLES_4 = dict(
    nodos_txt=widgets.Textarea(value="A, 0, 0\nB, 4, 0\nC, 2, 1.5", description="nodos (m):", continuous_update=False,
                               layout=_caja("110px")),
    barras_txt=widgets.Textarea(value="A, B\nA, C\nB, C", description="barras:", continuous_update=False,
                                layout=_caja("110px")),
    apoyos_txt=widgets.Textarea(value="A, x\nA, y\nB, y", description="apoyos:", continuous_update=False,
                                layout=_caja("90px")),
    cargas_txt=widgets.Textarea(value="C, 0, -10", description="cargas:", continuous_update=False, layout=_caja("70px")),
    unidad=widgets.Text(value="kN", description="unidad de fuerza:", continuous_update=False,
                        style={"description_width": "initial"}),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_armadura, **_CONTROLES_4);

# %%
# Casos de prueba de la sección 4 (resultado conocido)
_s5 = sp.sqrt(5)
_esperado_W = {"AB": sp.Rational(3, 2), "BC": sp.Rational(3, 2), "AD": -sp.Rational(3, 2) * _s5, "BD": sp.Rational(3, 2) * _s5,
               "BE": sp.Rational(3, 2) * _s5, "CE": -sp.Rational(3, 2) * _s5, "DE": -3}
_iguales = lambda a, b: set(a) == set(b) and all(sp.simplify(a[k] - b[k]) == 0 for k in a)
_cuadro = "A, 0, 0\nB, 2, 0\nC, 2, 2\nD, 0, 2"
PRUEBAS_4 = [
    ("triángulo del recorrido → AB = 20/3, AC = BC = -25/3, reacciones (Ax, Ay, By) = (0, 5, 5), determinada",
     lambda: _tri["fuerzas"] == {"AB": sp.Rational(20, 3), "AC": sp.Rational(-25, 3), "BC": sp.Rational(-25, 3)}
             and list(_tri["reacciones"].values()) == [0, 5, 5] and _tri["diagnostico"] == "determinada"
             and _tri["tipo"] == {"AB": "T", "AC": "C", "BC": "C"}
             and "punto decimal" in _mensaje(lambda: armadura("A, 0, 0\nB, 4, 0\nC, 2, 1,5", "A, B\nA, C\nB, C",
                                                              "A, x\nA, y\nB, y", "C, 0, -10"))),
    ("Warren → AB = BC = 3/2, AD = CE = -(3/2)√5, BD = BE = (3/2)√5, DE = -3, Ay = Cy = 3",
     lambda: _iguales(_warren["fuerzas"], _esperado_W) and _warren["reacciones"] == {"Ax": 0, "Ay": 3, "Cy": 3}),
    ("comprobacion_global: suma de fuerzas externas y momento respecto al primer nodo son cero en las pruebas 1 y 2",
     lambda: all(r["comprobacion_global"]["cero"] and r["comprobacion_global"]["suma_fuerzas"] == sp.zeros(2, 1)
                 and r["comprobacion_global"]["momento"] == 0 for r in (_tri, _warren))),
    ("triángulo A(0, 0), B(4, 0), C(2, 2) con tres apoyos verticales y carga en C → «inestable» (rango 5 < 6), sin traza",
     lambda: (lambda r: r["diagnostico"].startswith("inestable") and r["rango"] == 5 and r["conteo"] == (3, 3, 3)
                        and r["fuerzas"] is None)(
         armadura("A, 0, 0\nB, 4, 0\nC, 2, 2", "A, B\nA, C\nB, C", "A, y\nB, y\nC, y", "C, 0, -6"))
             and "INESTABLE" in _captura(calculadora_armadura, "A, 0, 0\nB, 4, 0\nC, 2, 2", "A, B\nA, C\nB, C",
                                         "A, y\nB, y\nC, y", "C, 0, -6", "kN", 4)[0]),
    ("cuadrado A(0, 0), B(2, 0), C(2, 2), D(0, 2) con sus 4 lados y apoyos A x, A y, B y → «mecanismo por conteo»",
     lambda: (lambda r: r["diagnostico"] == "mecanismo por conteo" and r["conteo"] == (4, 4, 3))(
         armadura(_cuadro, "A, B\nB, C\nC, D\nD, A", "A, x\nA, y\nB, y", "D, 1, 0"))),
    ("el mismo cuadrado con sus dos diagonales → «hiperestática por conteo»; con solo sus 4 lados y apoyos A x, A y, "
     "B x, B y, C y → «hiperestática por conteo e inestable (rango < 2n)»",
     lambda: (lambda r: r["diagnostico"] == "hiperestática por conteo" and r["conteo"] == (4, 6, 3))(
         armadura(_cuadro, "A, B\nB, C\nC, D\nD, A\nA, C\nB, D", "A, x\nA, y\nB, y", "D, 1, 0"))
             and (lambda r: r["diagnostico"] == "hiperestática por conteo e inestable (rango < 2n)"
                            and r["conteo"] == (4, 4, 5) and r["rango"] < 8 and r["fuerzas"] is None)(
                 armadura(_cuadro, "A, B\nB, C\nC, D\nD, A", "A, x\nA, y\nB, x\nB, y\nC, y", "D, 1, 0"))),
    ("linealidad: duplicar la carga de la Warren duplica todas las fuerzas y reacciones",
     lambda: (lambda r: _iguales(r["fuerzas"], {k: 2 * v for k, v in _esperado_W.items()})
                        and r["reacciones"] == {"Ax": 0, "Ay": 6, "Cy": 6})(
         armadura(_NODOS_W, _BARRAS_W, _APOYOS_W, "B, 0, -12"))),
    ("prueba de humo: la calculadora corre con sus valores iniciales y con la Warren, sin error y sin «Revisa la entrada»",
     lambda: _humo(calculadora_armadura, _CONTROLES_4)
             and _humo(calculadora_armadura, _CONTROLES_4, nodos_txt=_NODOS_W, barras_txt=_BARRAS_W,
                       apoyos_txt=_APOYOS_W, cargas_txt="B, 0, -6", unidad="N", n_cifras=6)),
]
corre_pruebas(PRUEBAS_4, 4)

# %% [markdown]
# **Contrasta (sección 4).** Para la armadura con $A(0,0)$ articulado, $B(2,0)$ sobre un rodillo y $C(1,1)$ con una carga de 4 kN hacia abajo, calcula a mano las fuerzas (empieza por el nodo $C$) y escribe `mis_fuerzas = (AB, AC, BC)` en kN, con tensión positiva.

# %%
mis_fuerzas = None     # tres números (AB, AC, BC) en kN, por ejemplo (1.5, -2.5, 0)

if mis_fuerzas is None:     # la referencia se muestra después de tu cálculo
    print("Falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda.")
else:
    try:
        mias = vector(mis_fuerzas, "mis_fuerzas")
    except ValueError as err:
        print(f"Revisa tu respuesta: {str(err).rstrip('.')}.")
    else:
        if mias.rows != 3:          # respuesta imposible: la referencia todavía no se muestra
            print(f"La armadura tiene tres barras, así que son tres fuerzas, AB, AC y BC (escribiste {mias.rows}): "
                  "revisa tu respuesta y vuelve a ejecutar la celda.")
        else:
            ref = armadura("A, 0, 0\nB, 2, 0\nC, 1, 1", "A, B\nA, C\nB, C", "A, x\nA, y\nB, y", "C, 0, -4")
            F = [ref["fuerzas"][k] for k in ("AB", "AC", "BC")]
            print("La calculadora da: " + ",  ".join(_magnitud(k, v, "kN", 4) for k, v in zip(("AB", "AC", "BC"), F)) + ".")
            if all(abs(float(a - b)) <= 0.01 for a, b in zip(mias, F)):
                print("Coinciden (con tolerancia de 0.01 kN).")
            elif all(abs(float(a + b)) <= 0.01 for a, b in zip(mias, F)):
                print("NO coinciden: tienes los valores correctos con los signos al revés. Revisa la convención: la "
                      "fuerza es positiva si la barra está en tensión (jala a sus nodos hacia ella).")
            else:
                print("NO coinciden: en C, las barras AC y BC forman 45° con la horizontal; por simetría son iguales, y "
                      "ΣFy = 0 da la suma de sus componentes verticales. Después, ΣFx = 0 en B da AB.")
            print("\nEl procedimiento completo:")
            _informe_armadura(ref, "kN", 4)
