# ID: ALG-U3-NB05
# Notebook: alg/u3_sistemas.ipynb · sección 3.5 interpretación geométrica en 2D y 3D
# Repositorio: alg/u3_sistemas/05_geometria.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 5. Interpretación geométrica en 2D y 3D
#
# Con dos incógnitas, cada ecuación $ax+by=c$ es una **recta** del plano; con tres, cada ecuación $ax+by+cz=d$ es un **plano** del espacio, con normal $\mathbf n=(a,b,c)$. La clasificación de la sección 4 tiene entonces una imagen, y los rangos dicen cuál.
#
# **Dos rectas** (sistema $2\times 2$):
#
# | $r$ | $r^*$ | Configuración |
# |:-:|:-:|---|
# | 2 | 2 | las rectas se cortan en un punto ($\det A\neq 0$) |
# | 1 | 1 | son la misma recta |
# | 1 | 2 | son paralelas y distintas |
#
# **Tres planos** (sistema $3\times 3$):
#
# | $r$ | $r^*$ | Configuración |
# |:-:|:-:|---|
# | 3 | 3 | un punto común |
# | 2 | 2 | una recta común: los tres planos pasan por una recta, o dos son el mismo y el tercero lo corta |
# | 1 | 1 | los tres son el mismo plano |
# | 2 | 3 | sin punto común: un **prisma** si ningún par de normales es paralelo; **dos planos paralelos cortados por un tercero** si exactamente un par lo es |
# | 1 | 2 | sin punto común: **planos paralelos** (al menos dos distintos) |
#
# Dos normales son paralelas si una es múltiplo de la otra, es decir, si $\mathbf n_i\times\mathbf n_j=\mathbf 0$. Con $r=1$ todas las normales son paralelas; con $r=2$ y $r^*=3$, la comparación de las normales separa el prisma de los planos paralelos cortados por un tercero.
#
# Con otra cantidad de ecuaciones vale la regla general de la sección 4: si $r<r^*$ no hay ningún punto común; si $r=r^*$, lo común es un punto ($r=n$), una recta ($n-r=1$) o un plano ($n-r=2$). Los planos de un sistema homogéneo pasan por el origen, así que siempre se cortan, al menos, en $\mathbf 0$.
#
# **Cómo leer la gráfica.** Cada ecuación tiene su propio estilo de línea (continua, discontinua, punteada...) además de su color, para que se distinga también en blanco y negro. El punto común es un círculo negro y la recta común, una línea negra gruesa. El título dice la configuración.

# %%
def _paralelas(u, v):
    """True si los vectores renglón u y v (no nulos) son paralelos: uno es múltiplo del otro."""
    return sp.Matrix.vstack(u, v).rank() == 1


def _pares(indices):
    return [(indices[a], indices[c]) for a in range(len(indices)) for c in range(a + 1, len(indices))]


def geometria(Ab, explicar=False, n_cifras=4):
    """Configuración geométrica de un sistema de 2 incógnitas (rectas) o 3 (planos), decidida con los rangos.

    Devuelve un diccionario con n, m, r, r_aum, clave (texto corto: 'punto', 'misma recta', 'paralelas', 'recta',
    'mismo plano', 'prisma', 'dos paralelos y un tercero', 'paralelos', 'sin punto común', 'todo' o 'sin solución'),
    configuracion (la frase completa), notas, el punto común y la recta común (particular, dirección) si existen."""
    Ab = aumentada(Ab)
    m, n = Ab.rows, Ab.cols - 1
    if n not in (2, 3):
        raise ValueError(f"la interpretación geométrica es para sistemas de 2 o 3 incógnitas (este tiene {_incognitas(n)}); "
                         "para clasificarlo usa la sección 4")
    A, b = Ab[:, :n], Ab[:, n]
    el = _elimina(Ab, n, con_b=True)
    r, r_aum = el["r"], el["r_aum"]
    cero = sp.zeros(1, n)
    nulas = [i for i in range(m) if A.row(i) == cero]               # ecuaciones sin incógnitas: 0 = c
    figuras = [i for i in range(m) if i not in nulas]               # las que sí son rectas o planos
    k = len(figuras)
    objeto, objetos = ("recta", "rectas") if n == 2 else ("plano", "planos")
    notas = []
    for i in nulas:
        notas.append(f"la ecuación {i + 1} dice 0 = {b[i]}: " + ("ningún punto la cumple" if b[i] != 0 else
                     f"no es {'una recta' if n == 2 else 'un plano'} y no restringe nada"))
    punto = recta = None
    if r == r_aum:
        _, xp, dirs, _ = solucion_general(el["R"], el["pivotes"], n)
        if r == n:
            punto = xp
        elif n - r == 1:
            recta = (xp, dirs[0])
    iguales = [(i, j) for i, j in _pares(figuras) if sp.Matrix.vstack(Ab.row(i), Ab.row(j)).rank() == 1]
    paralelas = [(i, j) for i, j in _pares(figuras) if _paralelas(A.row(i), A.row(j)) and (i, j) not in iguales]
    if any(b[i] != 0 for i in nulas):
        clave, texto = "sin solución", f"sin solución: una ecuación pide 0 = c con c ≠ 0 (r = {r}, r* = {r_aum})"
    elif k == 0:
        clave, texto = "todo", f"todas las ecuaciones son 0 = 0: cualquier punto {'del plano' if n == 2 else 'del espacio'} es solución"
    elif n == 2 and k == 2:
        if r == 2:
            clave, texto = "punto", f"las rectas se cortan en un punto, {_tupla(punto)}"
        elif r_aum == 1:
            clave, texto = "misma recta", "son la misma recta: todos sus puntos son solución"
        else:
            clave, texto = "paralelas", "son paralelas y distintas: no hay punto común"
    elif n == 3 and k == 3:
        if r == 3:
            clave, texto = "punto", f"un punto común: los tres planos se cortan en {_tupla(punto)}"
        elif r == r_aum == 2:
            clave = "recta"
            if iguales:
                i, j = iguales[0]
                texto = f"una recta común: los planos {i + 1} y {j + 1} son el mismo y el tercero lo corta en una recta"
            else:
                texto = "una recta común: los tres planos pasan por una misma recta"
        elif r == r_aum == 1:
            clave, texto = "mismo plano", "los tres son el mismo plano: todos sus puntos son solución"
        elif r == 2:
            if not paralelas:
                clave = "prisma"
                texto = ("prisma: no hay punto común; los planos se cortan de dos en dos en tres rectas paralelas "
                         "y ningún par de normales es paralelo")
            else:
                i, j = paralelas[0]
                clave = "dos paralelos y un tercero"
                texto = (f"dos planos paralelos cortados por un tercero: no hay punto común "
                         f"(los planos {i + 1} y {j + 1} tienen normales paralelas)")
        else:
            clave = "paralelos"
            texto = "planos paralelos: no hay punto común" + (" (dos de ellos son el mismo plano)" if iguales else
                                                              " (los tres son distintos)")
    elif k == 1:
        clave, texto = f"mism{'a' if n == 2 else 'o'} {objeto}", f"una sola ecuación: todos los puntos de su {objeto} son solución"
    elif n == 3 and k == 2:
        if r == 2:
            clave, texto = "recta", "una recta común: los dos planos se cortan en una recta"
        elif r_aum == 1:
            clave, texto = "mismo plano", "son el mismo plano: todos sus puntos son solución"
        else:
            clave, texto = "paralelos", "planos paralelos y distintos: no hay punto común"
    else:                                                       # más ecuaciones que las de la tabla: regla general
        if r < r_aum:
            if r == 1:
                clave, texto = (f"paralel{'as' if n == 2 else 'os'}",
                                f"{objetos} paralel{'as' if n == 2 else 'os'}, no {'todas iguales' if n == 2 else 'todos iguales'}: no hay punto común")
            else:
                clave, texto = "sin punto común", f"no hay ningún punto común a {'las' if n == 2 else 'los'} {k} {objetos}"
        elif r == n:
            clave, texto = "punto", f"{'las' if n == 2 else 'los'} {k} {objetos} pasan por un mismo punto, {_tupla(punto)}"
        elif r == 1:
            clave = f"mism{'a' if n == 2 else 'o'} {objeto}"
            texto = f"{'todas son la misma recta' if n == 2 else 'todos son el mismo plano'}: todos sus puntos son solución"
        else:
            clave, texto = "recta", f"una recta común: los {k} planos pasan por una misma recta"
    res = {"n": n, "m": m, "r": r, "r_aum": r_aum, "clave": clave, "configuracion": texto, "notas": notas,
           "punto": punto, "recta": recta, "paralelas": paralelas, "iguales": iguales}
    if explicar:
        _informe_geometria(Ab, res, n_cifras)
    return res


def _informe_geometria(Ab, res, n_cifras):
    _imprime_sistema(Ab)
    print(f"r = rango(A) = {res['r']},  r* = rango[A | b] = {res['r_aum']},  n = {res['n']}")
    print("Configuración:", res["configuracion"] + ".")
    for nota in res["notas"]:
        print("Nota:", nota + ".")
    if res["punto"] is not None and not all(e.is_integer for e in res["punto"]):
        print("Punto común ≈ (" + ", ".join(cifras(e, n_cifras) for e in res["punto"]) + ")")
    if res["recta"] is not None:
        p, d = res["recta"]
        print(f"Recta común: {_tupla(p)} + s·{_tupla(d)}, con s cualquier número (s es la variable libre).")
    if res["n"] == 3 and res["clave"] in ("prisma", "dos paralelos y un tercero"):
        A = Ab[:, :3]
        figuras = [i for i in range(Ab.rows) if A.row(i) != sp.zeros(1, 3)]
        for i, j in _pares(figuras):
            cruz = A.row(i).T.cross(A.row(j).T)
            print(f"  n{_sub(i + 1)} × n{_sub(j + 1)} = {_tupla(cruz)}" + ("  → normales paralelas" if cruz == sp.zeros(3, 1) else ""))


def _caja(Ab, res, minimo):
    """Centro y medio ancho de la ventana: alrededor del punto o la recta común, o de los puntos de cada figura
    más cercanos al origen."""
    n = res["n"]
    A, b = _flotante(Ab[:, :n]), _flotante(Ab[:, n]).ravel()
    cercanos = [b[i] / (A[i] @ A[i]) * A[i] for i in range(A.shape[0]) if A[i] @ A[i] > 0]
    if res["punto"] is not None:
        centro = _flotante(res["punto"]).ravel()
    elif res["recta"] is not None:
        p, d = (_flotante(v).ravel() for v in res["recta"])
        centro = p - (p @ d) / (d @ d) * d                       # el punto de la recta más cercano al origen
    elif cercanos:
        centro = np.mean(cercanos, axis=0)
    else:
        centro = np.zeros(n)
    lejos = max([float(np.linalg.norm(q - centro)) for q in cercanos] + [0.0])
    return centro, float(math.ceil(max(minimo, 1.4 * lejos + 1.5)))


_DE_CANTO = ("prisma", "dos paralelos y un tercero")


def _direccion_comun(Ab):
    """Dirección común exacta de los planos (n_i × n_j de dos normales que no son paralelas), o None si no hay."""
    A = Ab[:, :3]
    filas = [i for i in range(Ab.rows) if A.row(i) != sp.zeros(1, 3)]
    for i, j in _pares(filas):
        d = A.row(i).T.cross(A.row(j).T)
        if d != sp.zeros(3, 1):
            return d
    return None


def _vista(Ab, res):
    """(elevación, giro) automáticos para 3D. En un prisma, o con dos planos paralelos y un tercero, se mira a lo largo
    de la dirección común de los planos (calculada exacta): se ven de canto, como un triángulo o como dos rectas
    paralelas cortadas por otra. En los demás casos se busca, cerca de la vista usual (20°, -60°), una en la que ningún
    plano quede de canto."""
    A = _flotante(Ab[:, :3])
    normales = [f / np.linalg.norm(f) for f in A if np.linalg.norm(f) > 0]
    if not normales:
        return 20, -60
    if res["clave"] in _DE_CANTO:
        d = _direccion_comun(Ab)
        if d is None:
            return 20, -60
        d = _flotante(d).ravel()
        d = -d if d[2] < 0 else d
        return math.degrees(math.atan2(d[2], math.hypot(d[0], d[1]))), math.degrees(math.atan2(d[1], d[0]))
    vistas = []
    for elev in range(15, 45, 5):
        for azim in range(-180, 180, 5):
            if min(azim % 90, 90 - azim % 90) < 20:           # casi a lo largo de un eje: sus marcas se enciman
                continue
            e, a = math.radians(elev), math.radians(azim)
            ojo = np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)])
            de_canto = min(abs(ojo @ nn) for nn in normales)       # 0: algún plano se ve como una recta
            vistas.append((de_canto, abs(elev - 20) + min(abs(azim + 60), 360 - abs(azim + 60)), elev, azim))
    buenas = [v for v in vistas if v[0] >= 0.25]                  # ningún plano a menos de unos 15° de verse de canto
    elegida = min(buenas, key=lambda v: v[1]) if buenas else max(vistas)
    return elegida[2], elegida[3]


def _seccion(Ab):
    """Para un prisma o dos planos paralelos con un tercero: el corte con el plano d·x = 0, perpendicular a la
    dirección común d. Devuelve (centro, medio ancho, trazas) con los vértices de ese corte calculados exactos y la
    traza de cada plano (renglón, punto, dirección unitaria), que es la recta que se ve al mirar a lo largo de d."""
    A, b = Ab[:, :3], Ab[:, 3]
    d = _direccion_comun(Ab)
    filas = [i for i in range(Ab.rows) if A.row(i) != sp.zeros(1, 3)]
    vertices = []
    for i, j in _pares(filas):
        if A.row(i).T.cross(A.row(j).T) != sp.zeros(3, 1):        # dos planos que se cortan: su recta pasa por el corte
            M = sp.Matrix.vstack(A.row(i), A.row(j), d.T)
            vertices.append(_flotante(M.solve(sp.Matrix([b[i], b[j], 0]))).ravel())
    V = np.array(vertices)
    centro = V.mean(axis=0)
    L = max(1.0, 1.8 * float(np.max(np.linalg.norm(V - centro, axis=1))))
    dn = _flotante(d).ravel() / float(np.linalg.norm(_flotante(d)))
    trazas = []
    for i in filas:
        n_i = _flotante(A.row(i)).ravel()
        punto = centro + (float(b[i]) - n_i @ centro) / (n_i @ n_i) * n_i
        u = np.cross(dn, n_i)
        trazas.append((i, punto, u / np.linalg.norm(u)))
    return centro, L, trazas


def _dentro(p, u, bajo, alto):
    """Intervalo de t en el que p + t·u queda dentro de la caja [bajo, alto]."""
    t0, t1 = -np.inf, np.inf
    for k in range(len(p)):
        if abs(u[k]) > 1e-12:
            a, c = sorted(((bajo[k] - p[k]) / u[k], (alto[k] - p[k]) / u[k]))
            t0, t1 = max(t0, a), min(t1, c)
    return t0, t1


def graficar(Ab, giro=None, elevacion=None, n_cifras=4):
    """Dibuja las rectas (2 incógnitas) o los planos (3 incógnitas) del sistema, marca lo que tienen en común
    y escribe la configuración en el título. En 3D, sin giro ni elevación, la vista se elige sola (_vista); en un
    prisma o con dos planos paralelos y un tercero, cada plano se ve de canto y se dibuja como una recta gruesa
    con su estilo y su número."""
    Ab = aumentada(Ab)
    res = geometria(Ab)
    n, m = res["n"], Ab.rows
    A, b = _flotante(Ab[:, :n]), _flotante(Ab[:, n]).ravel()
    if not any(np.any(A[i] != 0) for i in range(m)):
        print("No hay nada que dibujar: todas las ecuaciones tienen sus coeficientes en cero, así que ninguna es "
              f"{'una recta' if n == 2 else 'un plano'}.")
        return res
    nombres = nombres_incognitas(n)
    centro, L = _caja(Ab, res, 4 if n == 2 else 3)
    titulo = f"r = {res['r']}, r* = {res['r_aum']}: {res['configuracion']}"
    if n == 2:
        fig, ax = plt.subplots(figsize=(6.2, 6.2))
        for i in range(m):
            a1, a2 = A[i]
            if a1 == 0 and a2 == 0:
                continue
            if abs(a2) >= abs(a1):
                xs = np.linspace(centro[0] - L, centro[0] + L, 2)
                ys = (b[i] - a1 * xs) / a2
            else:
                ys = np.linspace(centro[1] - L, centro[1] + L, 2)
                xs = (b[i] - a2 * ys) / a1
            ax.plot(xs, ys, linestyle=_ESTILOS[i], color=_COLORES[i], lw=2.2,
                    label=f"ecuación {i + 1}: {_ecuacion(Ab.row(i), nombres)}")
        if res["punto"] is not None:
            p = _flotante(res["punto"]).ravel()
            ax.plot(*p, "o", color="black", ms=9, zorder=5, label=f"punto común {_tupla(res['punto'])}")
        ax.set_xlim(centro[0] - L, centro[0] + L); ax.set_ylim(centro[1] - L, centro[1] + L); ax.set_aspect("equal")
        ax.axhline(0, color="black", lw=0.6); ax.axvline(0, color="black", lw=0.6)
        ax.grid(alpha=0.3); ax.set_xlabel("x"); ax.set_ylabel("y")
        lugar = "best"
    else:
        fig = plt.figure(figsize=(7.5, 6.8))
        ax = fig.add_subplot(projection="3d")
        de_canto = giro is None and elevacion is None and res["clave"] in _DE_CANTO and _direccion_comun(Ab) is not None
        if de_canto:                                             # cada plano se ve como una recta: una línea gruesa por plano
            centro, L, trazas = _seccion(Ab)
            bajo, alto = centro - L, centro + L
            for i, p, u in trazas:
                t0, t1 = _dentro(p, u, bajo, alto)
                P = p[None, :] + np.array([t0, t1])[:, None] * u[None, :]
                ax.plot(*P.T, linestyle=_ESTILOS[i], color=_COLORES[i], lw=3.5,
                        label=f"plano {i + 1}: {_ecuacion(Ab.row(i), nombres)}")
                ax.text(*(p + 0.92 * t1 * u), f" {i + 1}", fontsize=13, fontweight="bold", color="black")
            titulo += ". Vista a lo largo de la dirección común: cada plano se ve de canto, como una recta con su número"
        else:
            bajo, alto = centro - L, centro + L
            for i in range(m):
                normal = A[i]
                if not normal.any():
                    continue
                k = int(np.argmax(np.abs(normal)))             # se despeja la variable de mayor coeficiente
                u, v = [j for j in range(3) if j != k]
                U, V = np.meshgrid(np.linspace(bajo[u], alto[u], 9), np.linspace(bajo[v], alto[v], 9))
                W = (b[i] - normal[u] * U - normal[v] * V) / normal[k]
                W = np.where((W < bajo[k] - 1e-9) | (W > alto[k] + 1e-9), np.nan, W)   # fuera de la caja no se dibuja
                XYZ = [None, None, None]
                XYZ[u], XYZ[v], XYZ[k] = U, V, W
                ax.plot_surface(*XYZ, color=_COLORES[i], alpha=0.08, linewidth=0, shade=False)
                ax.plot_wireframe(*XYZ, color=_COLORES[i], linestyle=_ESTILOS[i], linewidth=1.2,
                                  label=f"plano {i + 1}: {_ecuacion(Ab.row(i), nombres)}")
            if res["punto"] is not None:
                p = _flotante(res["punto"]).ravel()
                ax.scatter(*p, s=90, color="black", depthshade=False, label=f"punto común {_tupla(res['punto'])}")
            if res["recta"] is not None:
                p, d = (_flotante(w).ravel() for w in res["recta"])
                t0, t1 = _dentro(p, d / np.linalg.norm(d), bajo, alto)
                if t0 < t1:
                    P = p[None, :] + np.array([t0, t1])[:, None] * (d / np.linalg.norm(d))[None, :]
                    ax.plot(*P.T, color="black", lw=3.5, label="recta común")
        ax.set_xlim(bajo[0], alto[0]); ax.set_ylim(bajo[1], alto[1]); ax.set_zlim(bajo[2], alto[2])
        ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
        ax.set_box_aspect((1, 1, 1))                         # misma escala en los tres ejes: los ángulos se ven bien
        ax.set_proj_type("ortho")                            # sin perspectiva, un plano de canto se ve como una recta
        auto = _vista(Ab, res)
        ax.view_init(elev=auto[0] if elevacion is None else elevacion, azim=auto[1] if giro is None else giro)
        lugar = "upper left"
    ax.set_title(textwrap.fill(titulo, 62), fontsize=10)
    if ax.get_legend_handles_labels()[0]:                    # solo si hay algo rotulado
        ax.legend(loc=lugar, fontsize=8)
    plt.show()
    return res

# %% [markdown]
# **Recorrido paso a paso.** Primero, tres planos con un solo punto común: $x+y+z=3$, $x-y+2z=2$, $2x+y-z=2$. Después, $x+y=2$, $y+z=2$, $x-z=1$: el primer renglón menos el segundo da $x-z=0$, que choca con la tercera ecuación ($x-z=1$), así que $r=2$ y $r^*=3$; ningún par de normales es paralelo, de modo que forman un prisma. La vista automática mira a lo largo de la dirección común de los tres planos, $\mathbf n_1\times\mathbf n_2$: los planos se ven de canto y cada uno se dibuja como una recta gruesa con su estilo y su número, y el prisma aparece como un triángulo. En la graficadora puedes desmarcar la vista automática y girarlo.

# %%
geometria("1, 1, 1, 3; 1, -1, 2, 2; 2, 1, -1, 2", explicar=True)
graficar("1, 1, 1, 3; 1, -1, 2, 2; 2, 1, -1, 2")
print()
geometria("1, 1, 0, 2; 0, 1, 1, 2; 1, 0, -1, 1", explicar=True)
graficar("1, 1, 0, 2; 0, 1, 1, 2; 1, 0, -1, 1");

# %% [markdown]
# **Graficadora.** Escribe un sistema de 2 incógnitas (rectas) o de 3 (planos) como matriz aumentada. Prueba también los de 3 incógnitas de los casos de prueba de abajo, por ejemplo `1, 0, 1, 2; 0, 1, 1, 3; 1, 1, 2, 5` (una recta común) o `1, 1, 1, 1; 1, 1, 1, 2; 1, 1, 1, 3` (planos paralelos). En 3D la vista se elige sola: en general, de modo que ningún plano quede de canto; en un prisma o con dos planos paralelos cortados por un tercero, a lo largo de la dirección común de los planos, que entonces se ven de canto, como rectas numeradas 1, 2, 3. Desmarca «vista automática» para girarla con los controles.

# %%
def calculadora_geometria(Ab_txt, vista_automatica, giro, elevacion, n_cifras):
    with bloque():                 # el texto aparece completo, de una sola vez; la gráfica va después
        try:
            Ab = aumentada(Ab_txt)
            res = geometria(Ab)
        except ValueError as err:
            print("Revisa la entrada:", err); return
        _informe_geometria(Ab, res, n_cifras)
    if vista_automatica:
        graficar(Ab)
    else:
        graficar(Ab, giro=giro, elevacion=elevacion)


widgets.interact(calculadora_geometria,
    Ab_txt=widgets.Textarea(value="2, 1, 5\n1, 3, 5", description="[A | b] =", continuous_update=False,
                            layout=widgets.Layout(width="340px", height="110px")),
    vista_automatica=widgets.Checkbox(value=True, description="vista automática (3D)"),
    giro=widgets.IntSlider(value=-60, min=-180, max=180, step=15, description="giro", continuous_update=False),
    elevacion=widgets.IntSlider(value=20, min=-10, max=90, step=5, description="elevación", continuous_update=False),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False));

# %%
# Casos de prueba de la sección 5 (resultado conocido)
PRUEBAS_5 = [
    ("2, 1, 5; 1, 3, 5 → las rectas se cortan en (2, 1)",
     lambda: (lambda g: g["clave"] == "punto" and g["punto"] == vector("2, 1")
                        and g["configuracion"].startswith("las rectas se cortan en un punto"))(geometria("2, 1, 5; 1, 3, 5"))),
    ("1, -2, 1; -2, 4, 3 → paralelas y distintas (r = 1, r* = 2)",
     lambda: (lambda g: g["clave"] == "paralelas" and (g["r"], g["r_aum"]) == (1, 2)
                        and "son paralelas y distintas" in g["configuracion"])(geometria("1, -2, 1; -2, 4, 3"))),
    ("1, -1, 1; 3, -3, 3 → la misma recta (r = r* = 1)",
     lambda: (lambda g: g["clave"] == "misma recta" and "son la misma recta" in g["configuracion"])(geometria("1, -1, 1; 3, -3, 3"))),
    ("1, 1, 1, 3; 1, -1, 2, 2; 2, 1, -1, 2 → un punto común, (1, 1, 1)",
     lambda: (lambda g: g["clave"] == "punto" and g["punto"] == vector("1, 1, 1")
                        and g["configuracion"].startswith("un punto común"))(geometria("1, 1, 1, 3; 1, -1, 2, 2; 2, 1, -1, 2"))),
    ("1, 0, 1, 2; 0, 1, 1, 3; 1, 1, 2, 5 → una recta común, que está en los tres planos",
     lambda: (lambda g: g["clave"] == "recta" and g["configuracion"].startswith("una recta común")
                        and all(matriz("1, 0, 1; 0, 1, 1; 1, 1, 2") * (g["recta"][0] + t * g["recta"][1]) == vector("2, 3, 5")
                                for t in (0, 1, -3)))(geometria("1, 0, 1, 2; 0, 1, 1, 3; 1, 1, 2, 5"))),
    ("1, 1, 0, 2; 0, 1, 1, 2; 1, 0, -1, 1 → prisma; y 1, 1, 1, 1; 1, 1, 1, 2; 1, 1, 1, 3 → planos paralelos",
     lambda: (lambda g1, g2: g1["clave"] == "prisma" and g1["configuracion"].startswith("prisma")
                             and g2["clave"] == "paralelos" and g2["configuracion"].startswith("planos paralelos"))(
         geometria("1, 1, 0, 2; 0, 1, 1, 2; 1, 0, -1, 1"), geometria("1, 1, 1, 1; 1, 1, 1, 2; 1, 1, 1, 3"))),
]
corre_pruebas(PRUEBAS_5, 5)

# %% [markdown]
# **Contrasta (sección 5).** Antes de graficar $x+y=2$, $2x+2y=5$, predice con los rangos qué vas a ver: calcula $r$ y $r^*$ y escribe una de estas tres frases: «se cortan en un punto», «son la misma recta» o «son paralelas y distintas».

# %%
mi_r = None              # tu rango de A, por ejemplo 2
mi_r_aum = None          # tu rango de [A | b], por ejemplo 1
mi_prediccion = None     # una de estas frases: "se cortan en un punto", "son la misma recta" o "son paralelas y distintas"

# Frases que se aceptan, ya normalizadas. Se compara la frase completa, así que una negación («no se cortan»)
# no cuenta como ninguna de las tres.
_PREDICCIONES = {
    "punto": ["se cortan en un punto", "se cortan", "se cortan en un solo punto", "un punto"],
    "misma recta": ["son la misma recta", "la misma recta", "misma recta"],
    "paralelas": ["son paralelas y distintas", "paralelas y distintas", "son paralelas", "paralelas"],
}


def _prediccion(texto):
    t = normaliza_texto(texto)
    for clave, frases in _PREDICCIONES.items():
        if t in frases:
            return clave
    raise ValueError(f"«{texto}» no es una de las tres frases: escribe «se cortan en un punto», «son la misma recta» "
                     "o «son paralelas y distintas»")


if mi_r is None or mi_r_aum is None or mi_prediccion is None:     # la referencia se muestra después de tu predicción
    print("falta tu predicción: escribe r, r* y lo que esperas ver, y vuelve a ejecutar la celda")
else:
    try:
        mio_r, mio_r_aum = numero(mi_r, "tu r"), numero(mi_r_aum, "tu r*")
        mia = _prediccion(mi_prediccion)
    except ValueError as err:
        print("Revisa tu respuesta:", err)
    else:
        if not all(v.is_integer and 0 <= v <= 2 for v in (mio_r, mio_r_aum)):   # respuesta imposible: sin referencia
            print("Con 2 ecuaciones y 2 incógnitas, r y r* son enteros de 0 a 2: revísalos y vuelve a ejecutar la celda.")
        else:
            ref = geometria("1, 1, 2; 2, 2, 5")
            print(f"La calculadora da: r = {ref['r']}, r* = {ref['r_aum']}; {ref['configuracion']}.")
            print("Rangos:", "coinciden" if (mio_r, mio_r_aum) == (ref["r"], ref["r_aum"]) else
                  "NO coinciden: escalona [A | b]; si el renglón 2 de A queda en ceros, mira qué queda en la columna de b.")
            print("Predicción:", "coincide" if mia == ref["clave"] else
                  "NO coincide: con r = 1 las rectas tienen la misma dirección; r* dice si además son la misma recta.")
            graficar("1, 1, 2; 2, 2, 5")
