# ID: ALG-U4-NB04
# Notebook: alg/u4_espacios.ipynb · sección 4.4 transformaciones lineales: rotación, escalamiento y proyección
# Repositorio: alg/u4_espacios/04_transformaciones.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 4. Transformaciones lineales: rotación, escalamiento y proyección
#
# **Transformación lineal.** Una función $T:\mathbb R^n\to\mathbb R^m$ es *lineal* si, para todos los vectores $\mathbf u,\mathbf v$ y todo número $c$,
#
# $$T(\mathbf u+\mathbf v)=T(\mathbf u)+T(\mathbf v)\qquad\text{y}\qquad T(c\,\mathbf u)=c\,T(\mathbf u).$$
#
# Entonces $T(\mathbf 0)=T(0\cdot\mathbf 0)=0\cdot T(\mathbf 0)=\mathbf 0$, y vale la **superposición**: $T(c_1\mathbf u_1+\dots+c_k\mathbf u_k)=c_1T(\mathbf u_1)+\dots+c_kT(\mathbf u_k)$.
#
# **Matriz estándar.** Como $\mathbf x=x_1\mathbf e_1+\dots+x_n\mathbf e_n$, la superposición da $T(\mathbf x)=x_1T(\mathbf e_1)+\dots+x_nT(\mathbf e_n)=A\mathbf x$ con
#
# $$A=[\,T(\mathbf e_1)\ \ T(\mathbf e_2)\ \cdots\ T(\mathbf e_n)\,]:$$
#
# basta saber a dónde van los vectores de la base estándar.
#
# **Las del plano y el espacio.**
#
# - **Rotación** un ángulo $\theta$ en sentido antihorario: $R(\theta)=\begin{bmatrix}\cos\theta&-\operatorname{sen}\theta\\ \operatorname{sen}\theta&\cos\theta\end{bmatrix}$. Ojo con la unidad: $\cos 60^\circ=\tfrac12$, pero $\cos(60\ \text{rad})\approx-0.952$. Por eso `rotacion(angulo, unidad)` te pide **siempre** la unidad, `"grados"` o `"radianes"`. En radianes puedes escribir el ángulo con pi: `"pi/6"`, `"2*pi/3"` o `"π/4"`.
# - **Rotaciones en $\mathbb R^3$** alrededor de un eje (antihorarias vistas desde la punta del eje, regla de la mano derecha):
#
#   $$R_x(\theta)=\begin{bmatrix}1&0&0\\0&\cos\theta&-\operatorname{sen}\theta\\0&\operatorname{sen}\theta&\cos\theta\end{bmatrix},\quad R_y(\theta)=\begin{bmatrix}\cos\theta&0&\operatorname{sen}\theta\\0&1&0\\-\operatorname{sen}\theta&0&\cos\theta\end{bmatrix},\quad R_z(\theta)=\begin{bmatrix}\cos\theta&-\operatorname{sen}\theta&0\\ \operatorname{sen}\theta&\cos\theta&0\\0&0&1\end{bmatrix}.$$
# - **Escalamiento:** $\operatorname{diag}(k_1,k_2)=\begin{bmatrix}k_1&0\\0&k_2\end{bmatrix}$ estira (o encoge) $k_1$ veces en $x$ y $k_2$ veces en $y$.
# - **Proyección** sobre la recta por el origen de dirección $\mathbf u\neq\mathbf 0$: $P=\dfrac{\mathbf u\,\mathbf u^T}{\mathbf u\cdot\mathbf u}$, porque $P\mathbf x=\dfrac{\mathbf u\cdot\mathbf x}{\mathbf u\cdot\mathbf u}\,\mathbf u$. Cumple $P^2=P$ (proyectar dos veces es proyectar una).
# - **Reflexión** respecto a esa recta: $H=2P-I$ (el punto medio entre $\mathbf x$ y $H\mathbf x$ es $P\mathbf x$). Cumple $H^2=I$. En el plano, $\det H=-1$; en $\mathbb R^n$, $\det H=(-1)^{n-1}$, porque la recta queda fija y las $n-1$ direcciones perpendiculares a ella cambian de sentido.
# - **Composición:** aplicar $T_1$ y después $T_2$ tiene matriz $A_2A_1$: **la que actúa primero va a la derecha**. En general $A_2A_1\neq A_1A_2$.
#
# **Núcleo e imagen.** $\operatorname{Nul}(T)=\{\mathbf x:\ T(\mathbf x)=\mathbf 0\}=\operatorname{Nul}(A)$ son los vectores que se aplastan en el origen; la imagen es $\operatorname{Col}(A)$, todo lo que $T$ alcanza. Por el teorema de la dimensión, $\dim\operatorname{Nul}(T)+\dim\operatorname{Im}(T)=n$. Una proyección del plano tiene núcleo de dimensión 1 (la recta perpendicular a $\mathbf u$) e imagen de dimensión 1 (la recta de $\mathbf u$); una rotación tiene núcleo $\{\mathbf 0\}$ e imagen $\mathbb R^2$.
#
# **No toda función "lineal" lo es.** Una traslación $\mathbf x\mapsto\mathbf x+\mathbf b$ con $\mathbf b\neq\mathbf 0$ y la función $y=mx+b$ de la prepa con $b\neq 0$ **no** son transformaciones lineales: no mandan $\mathbf 0$ a $\mathbf 0$. Por ejemplo, un transmisor de temperatura con salida de 4 a 20 mA, $I=4+0.16\,T$, no es lineal en este sentido (es *afín*), aunque su gráfica sea una recta.

# %%
_UNIDADES = {"grados": "grados", "grado": "grados", "radianes": "radianes", "radian": "radianes", "rad": "radianes"}


def _unidad(unidad):
    """'grados' o 'radianes' (acepta mayúsculas, acentos, el singular y 'rad'). No hay unidad por omisión: sin unidad,
    o con otra cosa, ValueError en español."""
    if unidad is None or (isinstance(unidad, str) and not unidad.strip()):
        raise ValueError("indica la unidad: «grados» o «radianes» (no hay unidad por omisión)")
    t = normaliza_texto(unidad) if isinstance(unidad, str) else ""
    if t not in _UNIDADES:
        raise ValueError(f"la unidad del ángulo debe ser «grados» o «radianes» (escribiste «{unidad}»)")
    return _UNIDADES[t]


_AYUDA_ANGULO = ("escribe un número (30, 37.5, 1/3) o, en radianes, un múltiplo de pi (pi/6, 2*pi/3, π/4, -pi)")


def angulo_exacto(angulo):
    """El ángulo como número exacto: 30 -> 30, '37.5' -> 75/2, 'pi/6' -> π/6, '2*pi/3', '2π/3', '-π/4'. El lector es
    restringido: solo números, pi, productos y cocientes (no evalúa expresiones libres)."""
    if not isinstance(angulo, str):
        return numero(angulo, "el ángulo")
    t = _normaliza_entrada(angulo).replace(" ", "").lower().replace("π", "pi")
    if "pi" not in t:
        return numero(angulo, "el ángulo")
    signo = 1
    if t[:1] in "+-":
        signo, t = (-1 if t[0] == "-" else 1), t[1:]
    t = re.sub(r"(?<=[\d.])pi", "*pi", t)                     # 2pi -> 2*pi
    trozos = re.split(r"([*/])", t)
    if not t or "" in trozos[::2]:
        raise ValueError(f"«{angulo}» no es un ángulo que entienda esta calculadora: {_AYUDA_ANGULO}")
    coef, potencia, op = sp.Integer(1), 0, "*"
    for i, trozo in enumerate(trozos):
        if i % 2:
            op = trozo
            continue
        if trozo == "pi":
            potencia += 1 if op == "*" else -1
            continue
        try:
            q = numero(trozo, "el ángulo")
        except ValueError:
            raise ValueError(f"«{angulo}» no es un ángulo que entienda esta calculadora: {_AYUDA_ANGULO}") from None
        if op == "/" and q == 0:
            raise ValueError(f"«{angulo}» divide entre cero")
        coef = coef * q if op == "*" else coef / q
    if potencia not in (0, 1):
        raise ValueError(f"«{angulo}» no es un ángulo que entienda esta calculadora: {_AYUDA_ANGULO}")
    return signo * coef * sp.pi ** potencia


def _angulo_rad(angulo, unidad):
    """El ángulo exacto en radianes: 30 grados -> pi/6; 30 radianes -> 30; 'pi/6' radianes -> pi/6."""
    u = _unidad(unidad)
    a = angulo_exacto(angulo)
    if u == "grados":
        if a.has(sp.pi):
            raise ValueError("con la unidad «grados» escribe el ángulo sin pi (por ejemplo 30); pi/6 ya está en radianes")
        return sp.rad(a)
    return a


def _num_angulo(a):
    """El número de un ángulo como texto: 30 -> '30', 75/2 -> '37.5', 1/3 -> '1/3', pi/6 -> 'π/6', 2*pi/3 -> '2π/3'."""
    if a.has(sp.pi):
        return _txt(a).replace("·π", "π")
    if a.is_integer:
        return str(a)
    q = int(a.q)
    while q % 2 == 0:
        q //= 2
    while q % 5 == 0:
        q //= 5
    return repr(float(a)) if q == 1 else str(a)          # decimal finito (37.5) o fracción (1/3)


def _texto_angulo(angulo, unidad):
    """'30 grados', '1 grado', 'π/6 radianes', '1 radián'."""
    a, u = angulo_exacto(angulo), _unidad(unidad)
    uno = abs(a) == 1
    nombre = ("grado" if uno else "grados") if u == "grados" else ("radián" if uno else "radianes")
    return f"{_num_angulo(a)} {nombre}"


def _tex_angulo(angulo, unidad):
    a = angulo_exacto(angulo)
    return rf"{_latex(a)}^\circ" if _unidad(unidad) == "grados" else rf"{_latex(a)}\ \text{{rad}}"


def rotacion(angulo, unidad=None):
    """R(θ) exacta, antihoraria. La unidad es obligatoria: rotacion(30, "grados"), rotacion("pi/6", "radianes");
    sin unidad, ValueError «indica la unidad»."""
    t = _angulo_rad(angulo, unidad)
    c, s = sp.cos(t), sp.sin(t)
    return sp.Matrix([[c, -s], [s, c]])


def rotacion3(eje, angulo, unidad=None):
    """Rotación exacta de R³ alrededor del eje 'x', 'y' o 'z' (antihoraria vista desde la punta del eje, regla de la
    mano derecha). La unidad es obligatoria, como en rotacion."""
    e = normaliza_texto(eje) if isinstance(eje, str) else ""
    if e not in ("x", "y", "z"):
        raise ValueError(f"el eje de rotación debe ser «x», «y» o «z» (escribiste «{eje}»)")
    t = _angulo_rad(angulo, unidad)
    c, s = sp.cos(t), sp.sin(t)
    if e == "x":
        return sp.Matrix([[1, 0, 0], [0, c, -s], [0, s, c]])
    if e == "y":
        return sp.Matrix([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    return sp.Matrix([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def escalamiento(*k):
    """diag(k1, k2, ...) exacta: escalamiento(1.5, 0.5) -> diag(3/2, 1/2)."""
    if not k:
        raise ValueError("escribe al menos un factor de escala, por ejemplo escalamiento(2, 1)")
    if len(k) > 6:
        raise ValueError(f"escribiste {len(k)} factores; esta calculadora acepta hasta 6")
    return sp.diag(*[numero(v, f"el factor k{i}") for i, v in enumerate(k, 1)])


def _direccion(u):
    """La dirección u de una recta: un vector no nulo de hasta 6 componentes."""
    u = vector(u, "la dirección u")
    if u.rows > 6:
        raise ValueError(f"la dirección u tiene {u.rows} componentes; esta calculadora acepta hasta 6")
    if u == sp.zeros(u.rows, 1):
        raise ValueError("la dirección u no puede ser el vector cero: el 0 no define ninguna recta")
    return u


def proyeccion(u):
    """Proyección exacta sobre la recta por el origen de dirección u: P = u uᵀ / (u·u)."""
    u = _direccion(u)
    return u * u.T / (u.T * u)[0]


def reflexion(u):
    """Reflexión exacta respecto a la recta por el origen de dirección u: H = 2P - I."""
    P = proyeccion(u)
    return 2 * P - sp.eye(P.rows)


def matriz_desde_imagenes(imagenes):
    """Matriz estándar A = [T(e1) ... T(en)] a partir de las imágenes de la base estándar, escritas una por renglón
    (como lista de vectores): '2, 1; -1, 3' son T(e1) = (2, 1) y T(e2) = (-1, 3), y quedan como columnas."""
    return lista_vectores(imagenes, "la lista de imágenes T(e1), ..., T(en)")


def _matriz_t(A, nombre="A"):
    """Matriz de una transformación: texto POR RENGLONES (racional, hasta 6×6) o Matrix exacta de sympy, que puede
    tener raíces, senos y cosenos (la de una rotación)."""
    if isinstance(A, sp.MatrixBase) and not all(sp.sympify(e).is_Rational for e in A):
        M = sp.Matrix(A)
        _revisa_tamano(M, nombre)
        if not all(e.is_number and e.is_real for e in M):
            raise ValueError(f"{nombre} debe tener entradas numéricas reales")
        return M
    return matriz_limitada(A, nombre)


def nucleo_imagen(A):
    """Núcleo e imagen de T(x) = A x: diccionario con rango, base_nucleo (base exacta de Nul(A)), base_imagen
    (columnas pivote de la A ORIGINAL), dim_nucleo, dim_imagen, n (columnas: T va de Rⁿ...) y m (renglones: ...a Rᵐ)."""
    A = _matriz_t(A)
    m, n = A.shape
    el = _elimina(A, n, con_b=False)
    piv = el["pivotes"]
    nucleo = _direcciones(el["R"], piv, n)
    return {"A": A, "rango": el["r"], "base_nucleo": nucleo, "base_imagen": [A[:, c] for c in piv],
            "dim_nucleo": len(nucleo), "dim_imagen": el["r"], "n": n, "m": m, "pivote": [c + 1 for c in piv]}


def _iguales(a, b):
    """True si los vectores a y b son iguales: exacto con fracciones, raíces, senos...; con decimales de punto
    flotante (si la función los usa), con una tolerancia relativa de 1e-9."""
    if a.shape != b.shape:
        return False
    for p, q in zip(a, b):
        d = sp.sympify(p - q)
        if d.has(sp.Float):
            if abs(float(d)) > 1e-9 * max(1.0, abs(float(p)), abs(float(q))):
                return False
        elif not _es_cero(d):
            return False
    return True


_SIN_CONTRAEJEMPLO = ("no se encontró contraejemplo con estas pruebas (es una prueba con vectores fijos, no una "
                      "demostración); si T es lineal, su matriz estándar es A = [T(e1) ... T(en)]")


def prueba_linealidad(f, n):
    """Busca un contraejemplo a la linealidad de f: Rⁿ → Rᵐ. f recibe una Matrix columna de sympy de tamaño n y
    devuelve otra (o una lista). Revisa primero T(0) = 0; después la aditividad T(u + v) = T(u) + T(v) y la
    homogeneidad T(c·u) = c·T(u) con vectores fijos, enteros y no enteros (e1, ..., en, (1, 2, ..., n),
    (3, -1, 4, ...), (1/2, -1/3, 2/5, ...)) y los escalares 2, -3, 1/2 y -1/3, en forma exacta; al final, que
    T(u) = A·u con A = [T(e1) ... T(en)]. Devuelve un diccionario con lineal (False si encontró un contraejemplo),
    falla (el primer contraejemplo, como texto), mensaje y, si no encontró contraejemplo, matriz. Es una PRUEBA, no
    una demostración: puede descubrir que f no es lineal, pero no garantizar que lo sea."""
    if isinstance(n, bool) or not isinstance(n, (int, np.integer)) or not 1 <= int(n) <= 6:
        raise ValueError("n, el número de componentes del vector de entrada, debe ser un entero de 1 a 6")
    n = int(n)
    tam = [None]

    def T(v):
        try:
            y = f(v.copy())
            y = sp.Matrix(y) if isinstance(y, (list, tuple, sp.MatrixBase, np.ndarray)) else sp.Matrix([y])
        except Exception as err:
            raise ValueError(f"la función no se pudo evaluar en x = {_tupla(v)} ({type(err).__name__}: {err})") from None
        y = y.reshape(len(y), 1)
        if len(y) == 0 or not all(sp.sympify(e).is_number for e in y):
            raise ValueError(f"T{_tupla(v)} = {_tupla(y)} no es un vector de números")
        if tam[0] is None:
            tam[0] = y.rows
        elif y.rows != tam[0]:
            raise ValueError(f"la función da vectores de tamaños distintos ({tam[0]} y {y.rows} componentes)")
        return y

    def falla(texto):
        return {"lineal": False, "falla": texto, "mensaje": f"no es lineal: {texto}", "matriz": None}

    T0 = T(sp.zeros(n, 1))
    if not _iguales(T0, sp.zeros(T0.rows, 1)):
        return falla(f"T(0) = {_tupla(T0)}, no es 0")
    base = [sp.eye(n)[:, j] for j in range(n)]
    prueba = list(base)
    R = sp.Rational
    for v in (sp.Matrix(list(range(1, n + 1))), sp.Matrix([3, -1, 4, -2, 5, -3][:n]),
              sp.Matrix([R(1, 2), R(-1, 3), R(2, 5), R(-1, 7), R(3, 11), R(-2, 13)][:n])):
        if v not in prueba:
            prueba.append(v)
    for i, u in enumerate(prueba):
        for w in prueba[i + 1:]:
            izq, der = T(u + w), T(u) + T(w)
            if not _iguales(izq, der):
                return falla(f"T(u + v) ≠ T(u) + T(v) con u = {_tupla(u)} y v = {_tupla(w)}: "
                             f"T(u + v) = {_tupla(izq)} y T(u) + T(v) = {_tupla(der)}")
    for c in (2, -3, R(1, 2), R(-1, 3)):
        for u in prueba:
            izq, der = T(c * u), c * T(u)
            if not _iguales(izq, der):
                cc = _paren_txt(c)
                return falla(f"T({cc}·u) ≠ {cc}·T(u) con u = {_tupla(u)}: T({cc}·u) = {_tupla(izq)} y "
                             f"{cc}·T(u) = {_tupla(der)}")
    A = sp.Matrix.hstack(*[T(e) for e in base])
    for u in prueba:
        if not _iguales(T(u), A * u):
            return falla(f"T(u) ≠ A·u con u = {_tupla(u)}, A = [T(e1) ... T(en)]")
    return {"lineal": True, "falla": None, "mensaje": _SIN_CONTRAEJEMPLO, "matriz": A}


def _paren_txt(c):
    """Un escalar para escribir delante de '·u': 2 -> '2', 1/2 -> '(1/2)', -3 -> '(-3)'."""
    return str(c) if (c.is_integer and c >= 0) else f"({c})"


# --- dibujo de una transformación del plano -----------------------------------------------------
_CASA = np.array([(0, 0), (2, 0), (2, 1.5), (1, 2.3), (0, 1.5)], dtype=float)
_CUADRADO = np.array([(0, 0), (1, 0), (1, 1), (0, 1)], dtype=float)


def _dibuja_en(ax, A, figura="casa", titulo=None):
    """Dibuja en ax la figura original (punteada) y su imagen bajo A (rellena gris), e1 y T(e1) (rotulados en la
    leyenda, para que el rótulo no se encime en la figura); si rango < 2, segmentos finos discontinuos en la
    dirección del núcleo: los puntos de cada uno caen en el mismo lugar."""
    Af = _flotante(A)
    F = _CASA if figura == "casa" else _CUADRADO
    G = F @ Af.T
    ni = nucleo_imagen(A)
    cerrada = lambda P: np.vstack([P, P[:1]])
    ax.fill(*G.T, facecolor="0.8", edgecolor="black", lw=1.8, hatch="..", alpha=0.9, label=f"imagen de la {figura}")
    ax.plot(*cerrada(G).T, "-", color="black", lw=1.8)
    ax.plot(*cerrada(F).T, ":", color="black", lw=2.0, label=f"{figura} original")
    e1, Te1 = np.array([1.0, 0.0]), Af @ np.array([1.0, 0.0])
    L = 1.2 * max(1.0, float(np.abs(np.vstack([F, G, Te1])).max())) + 0.4
    if ni["rango"] == 1:
        d = _flotante(ni["base_nucleo"][0]).ravel()
        d = d / np.linalg.norm(d)
        for k, p in enumerate(F):
            ax.plot(*np.array([p - 3 * L * d, p + 3 * L * d]).T, "--", color="0.45", lw=0.8,
                    label="dirección del núcleo: cada línea cae en un solo punto" if k == 0 else None)
    _flecha(ax, (0, 0), e1, _COLORES[0], "--", 2.2)
    _flecha(ax, (0, 0), Te1, _COLORES[1], "-", 2.6)
    ax.plot([], [], "--", color=_COLORES[0], lw=2.2, label="e1 (flecha discontinua)")
    ax.plot([], [], "-", color=_COLORES[1], lw=2.6, label="T(e1) (flecha continua)")
    ax.set_xlim(-L, L); ax.set_ylim(-L, L); ax.set_aspect("equal")
    ax.axhline(0, color="black", lw=0.6); ax.axvline(0, color="black", lw=0.6)
    ax.grid(alpha=0.4); ax.set_xlabel("x"); ax.set_ylabel("y")
    if titulo is None:
        titulo = {2: "rango 2: la figura cambia de tamaño, forma u orientación, pero no se aplasta",
                  1: "rango 1: la figura se aplasta sobre una recta (el núcleo es una recta)",
                  0: "rango 0: todo va a dar al origen"}[ni["rango"]]
    ax.set_title(textwrap.fill(titulo, 55), fontsize=10)
    ax.legend(loc="lower left", fontsize=8)


def dibuja_transformacion(A, figura="casa"):
    """Dibuja en el plano una casa (o un cuadrado) y su imagen bajo la transformación de matriz A (de 2×2)."""
    A = _matriz_t(A)
    if A.shape != (2, 2):
        raise ValueError(f"el dibujo es para transformaciones del plano: A debe ser de 2×2 (esta es {dims(A)})")
    if figura not in ("casa", "cuadrado"):
        raise ValueError("la figura debe ser «casa» o «cuadrado»")
    fig, ax = plt.subplots(figsize=(6.2, 6.2))
    _dibuja_en(ax, A, figura)
    plt.show()

# %% [markdown]
# **Recorrido paso a paso.** Primero, $R(30^\circ)$ aplicada a $(2,0)$: el resultado exacto es $(\sqrt3,1)$. Después, la proyección sobre la recta de $\mathbf u=(3,4)$: $P=\tfrac1{25}\begin{bmatrix}9&12\\12&16\end{bmatrix}$; $(5,0)$ va a $(9/5,12/5)$, $(4,-3)$ (perpendicular a $\mathbf u$) va al origen y $P^2=P$. La reflexión respecto a esa recta manda $(5,0)$ a $(-7/5,24/5)$. Al final, dos rotaciones seguidas: $R(20^\circ)R(25^\circ)=R(45^\circ)$, y los dibujos de la casa bajo la rotación y bajo la proyección (que la aplasta sobre una recta).

# %%
_R30 = rotacion(30, "grados")
muestra("R(30 grados)", _R30, 4, tex=r"R(30^\circ)")
print("R(30 grados)·(2, 0) =", _tupla(_R30 * vector("2, 0"), 4))
print()
_P = proyeccion("3, 4")
muestra("P", _P, 4)
print("P = (1/25)·" + _en_linea(25 * _P))
print("P·(5, 0) =", _tupla(_P * vector("5, 0"), 4), "  P·(4, -3) =", _tupla(_P * vector("4, -3")),
      "  ¿P² = P?", "sí" if _P * _P == _P else "no")
_H = reflexion("3, 4")
muestra("H = 2P - I", _H, 4)
print("H·(5, 0) =", _tupla(_H * vector("5, 0"), 4), "  det H =", _H.det(),
      "  ¿H² = I?", "sí" if _H * _H == sp.eye(2) else "no")
print()
_C = sp.simplify(rotacion(20, "grados") * rotacion(25, "grados"))
muestra("R(20 grados)·R(25 grados)", _C, 4, tex=r"R(20^\circ)\,R(25^\circ)")
print("¿R(20 grados)·R(25 grados) = R(45 grados)?", "sí" if _C == rotacion(45, "grados") else "no")
dibuja_transformacion(_R30)
dibuja_transformacion(_P);

# %% [markdown]
# **¿Es lineal?** `prueba_linealidad(f, n)` busca un contraejemplo con vectores enteros y no enteros. La función del transmisor de 4 a 20 mA, $I=4+0.16\,T$, falla desde $T(\mathbf 0)$; $T(x,y)=(x^2,y)$ falla la aditividad; $T(x)=\lfloor x\rfloor$ (parte entera) pasa las pruebas con enteros, pero no con $\tfrac12$. $T(x,y)=(2x-y,\ x+3y)$ no da ningún contraejemplo y la función devuelve su matriz estándar; ojo: es una **prueba**, no una demostración (la demostración es verificar las dos propiedades con letras). La misma matriz sale de `matriz_desde_imagenes`, con $T(\mathbf e_1)$ y $T(\mathbf e_2)$ escritos uno por renglón (como lista de vectores: cada renglón se vuelve una columna).

# %%
print(prueba_linealidad(lambda v: sp.Matrix([4 + sp.Rational(4, 25) * v[0]]), 1)["mensaje"])
print(prueba_linealidad(lambda v: sp.Matrix([v[0] ** 2, v[1]]), 2)["mensaje"])
print(prueba_linealidad(lambda v: sp.Matrix([sp.floor(v[0])]), 1)["mensaje"])
_lin = prueba_linealidad(lambda v: sp.Matrix([2 * v[0] - v[1], v[0] + 3 * v[1]]), 2)
print(_lin["mensaje"][0].upper() + _lin["mensaje"][1:] + ":", _en_linea(_lin["matriz"]))
_A_img = matriz_desde_imagenes("2, 1; -1, 3")          # T(e1) = (2, 1) y T(e2) = (-1, 3), uno por renglón
print("Con las imágenes T(e1) = (2, 1) y T(e2) = (-1, 3) como columnas:", _en_linea(_A_img),
      "(la misma)" if _A_img == _lin["matriz"] else "(¡distinta!)")

# %% [markdown]
# **Calculadora.** Elige el tipo de transformación y sus datos: para la rotación, el ángulo **y su unidad** (en radianes puedes escribir `pi/6`, `2*pi/3` o `π/4`); para el escalamiento, $k_1$ y $k_2$; para la proyección y la reflexión, la dirección $\mathbf u$ de la recta; para la matriz propia, la matriz **por renglones**. Los controles que no usa el tipo elegido se ignoran. Prueba la rotación con 30 en radianes y compara con 30 grados.

# %%
_TIPOS = ["rotación", "escalamiento", "proyección sobre una recta", "reflexión respecto a una recta", "matriz propia"]


def _tex_tupla(v):
    return r"\left(" + ",\\ ".join(_latex(e) for e in v) + r"\right)"


def _transformacion(tipo, angulo, unidad, k1, k2, u_txt, M_txt):
    """(A, nombre, tex, avisos) de la transformación elegida en la calculadora."""
    avisos = []
    if tipo == "rotación":
        u_ang = _unidad(unidad)
        A = rotacion(angulo, u_ang)
        nombre, tex = f"R({_texto_angulo(angulo, u_ang)})", rf"R\left({_tex_angulo(angulo, u_ang)}\right)"
        if u_ang == "radianes" and abs(angulo_exacto(angulo)) > 2 * sp.pi:
            avisos.append(f"escribiste {_texto_angulo(angulo, u_ang)}, más de una vuelta completa (2π ≈ 6.283 radianes). "
                          "¿Seguro que está en radianes? Si lo pensaste en grados, cambia la unidad a «grados».")
    elif tipo == "escalamiento":
        A = escalamiento(k1, k2)
        nombre = f"S = diag({_txt(A[0, 0])}, {_txt(A[1, 1])})"
        tex = rf"S=\operatorname{{diag}}\left({_latex(A[0, 0])},\ {_latex(A[1, 1])}\right)"
    elif tipo in ("proyección sobre una recta", "reflexión respecto a una recta"):
        u = _direccion(u_txt)
        if u.rows != 2:
            raise ValueError(f"en el plano, la dirección u debe tener 2 componentes (escribiste {u.rows})")
        if tipo.startswith("proyección"):
            A, letra, frase = proyeccion(u), "P", "proyección sobre la recta de"
        else:
            A, letra, frase = reflexion(u), "H", "reflexión respecto a la recta de"
        nombre = f"{letra} ({frase} u = {_tupla(u)})"
        tex = rf"{letra}\ \text{{({frase} }}\mathbf u={_tex_tupla(u)}\text{{)}}"
    elif tipo == "matriz propia":
        A, nombre, tex = matriz_limitada(M_txt, "la matriz propia"), "A (matriz propia)", r"A\ \text{(matriz propia)}"
    else:
        raise ValueError(f"tipo desconocido: «{tipo}»; elige uno de: {', '.join(_TIPOS)}")
    return A, nombre, tex, avisos


def _vec_texto(v, n_cifras):
    """Un vector exacto con sus decimales; si lo exacto es muy largo (raíces anidadas), solo los decimales."""
    return _tupla(v, n_cifras) if _exacta_legible(v) else "≈ " + _aprox(v, n_cifras)


def _vec_base(v, n_cifras):
    """Un vector de una base, para mostrarlo: si es racional, escalado a enteros sin factor común (como en las
    secciones 1 y 2); si no, exacto o en decimales."""
    if all(e.is_Rational for e in v) and any(e != 0 for e in v):
        return _tupla(_primitivo(v))
    return _vec_texto(_simplifica(v), n_cifras)


def _informe_transformacion(A, nombre, x, n_cifras, tex=None):
    m, n = A.shape
    corto = nombre.split(" (")[0].split(" = ")[0]
    muestra(nombre, A, n_cifras, tex=tex)
    print(f"T va de {_R(n)} a {_R(m)}: T(x) = {corto}·x.")
    if m == n:
        d = sp.simplify(A.det())
        if not _exacta_legible(sp.Matrix([d])):
            print(f"det {corto} ≈ {cifras(d, n_cifras)}")
        else:
            print(f"det {corto} = {_txt(d)}" + ("" if d.is_Rational else f" ≈ {cifras(d, n_cifras)}"))
    else:
        print("det: solo existe para matrices cuadradas.")
    ni = nucleo_imagen(A)
    print(f"rango = {ni['rango']}")
    nucleo = ni["base_nucleo"]
    if not nucleo:
        print("Núcleo: {0} (dim 0): solo el 0 va al 0.")
    else:
        print(f"Núcleo: gen{{{', '.join(_vec_base(v, n_cifras) for v in nucleo)}}} (dim {len(nucleo)}): " +
              ("ese vector y sus múltiplos van al 0." if len(nucleo) == 1 else
               "esos vectores y sus combinaciones van al 0."))
    if ni["dim_imagen"] == m:
        print(f"Imagen: todo {_R(m)} (dim {m}).")
    elif not ni["base_imagen"]:
        print("Imagen: {0} (dim 0).")
    else:
        print(f"Imagen: gen{{{', '.join(_vec_base(v, n_cifras) for v in ni['base_imagen'])}}} (dim {ni['dim_imagen']}).")
    print(f"Teorema de la dimensión: dim Nul + dim Im = {ni['dim_nucleo']} + {ni['dim_imagen']} = {n} = n")
    Tx = _simplifica(A * x)
    print(f"T(x) = {corto}·{_tupla(x)} " + ("= " if _exacta_legible(Tx) else "") + _vec_texto(Tx, n_cifras))
    return ni


def calculadora_transformacion(tipo, angulo, unidad, k1, k2, u_txt, M_txt, x_txt, n_cifras):
    with bloque():                 # el texto aparece completo, de una sola vez; el dibujo va después
        try:
            A, nombre, tex, avisos = _transformacion(tipo, angulo, unidad, k1, k2, u_txt, M_txt)
            x = vector(x_txt, "x")
            if x.rows != A.cols:
                raise ValueError(f"x tiene {_componentes(x.rows)} y la matriz tiene {A.cols} "
                                 f"{'columna' if A.cols == 1 else 'columnas'}: deben coincidir")
        except ValueError as err:
            print("Revisa la entrada:", err); return
        for aviso in avisos:
            print("Aviso:", aviso)
        _informe_transformacion(A, nombre, x, n_cifras, tex)
    if A.shape == (2, 2):
        dibuja_transformacion(A)
    else:
        print(f"El dibujo de la casa es para matrices de 2×2 (transformaciones del plano); esta es {dims(A)}.")


_CONTROLES_4A = dict(
    tipo=widgets.Dropdown(options=_TIPOS, value="rotación", description="tipo:"),
    angulo=widgets.Text(value="30", description="ángulo:", continuous_update=False),
    unidad=widgets.ToggleButtons(options=["grados", "radianes"], value="grados", description="unidad:"),
    k1=widgets.FloatText(value=1.5, description="k1:", continuous_update=False),
    k2=widgets.FloatText(value=0.5, description="k2:", continuous_update=False),
    u_txt=widgets.Text(value="3, 4", description="u =", continuous_update=False),
    M_txt=widgets.Textarea(value="1, 1; 0, 1", description="matriz:", continuous_update=False,
                           layout=widgets.Layout(width="300px", height="80px")),
    x_txt=widgets.Text(value="2, 0", description="x =", continuous_update=False),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_transformacion, **_CONTROLES_4A);

# %% [markdown]
# **Calculadora de composición.** Elige $T_1$ (actúa primero) y $T_2$ (actúa después). La matriz de «primero $T_1$, después $T_2$» es $A_2A_1$; la calculadora la compara con $A_1A_2$ y dibuja la casa con los dos órdenes.

# %%
_SENCILLAS = {
    "rotación de 45°": lambda: rotacion(45, "grados"),
    "rotación de -30°": lambda: rotacion(-30, "grados"),
    "escalamiento diag(2, 1)": lambda: escalamiento(2, 1),
    "corte [[1, 1], [0, 1]]": lambda: matriz("1, 1; 0, 1"),
    "reflexión respecto al eje y": lambda: reflexion("0, 1"),
    "proyección sobre el eje y": lambda: proyeccion("0, 1"),
}


def calculadora_composicion(T1, T2, n_cifras):
    A1, A2 = _SENCILLAS[T1](), _SENCILLAS[T2]()
    P21, P12 = _simplifica(A2 * A1), _simplifica(A1 * A2)
    igual = all(_es_cero(e) for e in P21 - P12)
    with bloque():
        print(f"T1 = {T1} (actúa primero), T2 = {T2} (actúa después).")
        muestra("A1", A1, n_cifras)
        muestra("A2", A2, n_cifras)
        muestra("A2·A1 (primero T1, después T2)", P21, n_cifras,
                tex=r"A_2A_1\ \text{(primero }T_1\text{, después }T_2\text{)}")
        muestra("A1·A2 (primero T2, después T1)", P12, n_cifras,
                tex=r"A_1A_2\ \text{(primero }T_2\text{, después }T_1\text{)}")
        print("Coinciden: con estas dos, el orden no importa." if igual else
              "NO coinciden: el orden importa (el producto de matrices no es conmutativo).")
    fig, ejes = plt.subplots(1, 2, figsize=(11, 5.6))
    _dibuja_en(ejes[0], P21, titulo="primero T1, después T2: A2·A1")
    _dibuja_en(ejes[1], P12, titulo="primero T2, después T1: A1·A2")
    plt.tight_layout()
    plt.show()


_CONTROLES_4B = dict(
    T1=widgets.Dropdown(options=list(_SENCILLAS), value="rotación de 45°", description="T1 (primero):",
                        style={"description_width": "initial"}),
    T2=widgets.Dropdown(options=list(_SENCILLAS), value="reflexión respecto al eje y", description="T2 (después):",
                        style={"description_width": "initial"}),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_composicion, **_CONTROLES_4B);

# %%
# Casos de prueba de la sección 4 (resultado conocido)
_x20 = vector("2, 0")
_NO_LINEALES = [(lambda v: sp.Matrix([sp.floor(v[0])]), 1),                       # parte entera
                (lambda v: sp.Matrix([v[0] + sp.sin(2 * sp.pi * v[0])]), 1),       # x + sen(2πx)
                (lambda v: sp.Matrix([v[0] * sp.cos(sp.pi * v[1]) ** 2]), 2)]      # x·cos²(πy)
PRUEBAS_4 = [
    ("R(90 grados) = [[0, -1], [1, 0]]",
     lambda: rotacion(90, "grados") == sp.Matrix([[0, -1], [1, 0]])),
    ("R(30 grados)·(2, 0) = (√3, 1) exacto",
     lambda: rotacion(30, "grados") * _x20 == sp.Matrix([sp.sqrt(3), 1])),
    ("R(30 radianes) ≠ R(30 grados) y su entrada [0, 0] ≈ 0.1543 (cos 30 rad)",
     lambda: rotacion(30, "radianes") != rotacion(30, "grados")
             and cifras(rotacion(30, "radianes")[0, 0], 4) == "0.1543"),
    ("R(25 grados)·R(35 grados) = R(60 grados)",
     lambda: sp.simplify(rotacion(25, "grados") * rotacion(35, "grados") - rotacion(60, "grados")) == sp.zeros(2, 2)
             and np.allclose(_flotante(rotacion(25, "grados") * rotacion(35, "grados")),
                             _flotante(rotacion(60, "grados")), atol=1e-12, rtol=0)),
    ("proyección sobre u = (3, 4): P = (1/25)[[9, 12], [12, 16]]; P(5, 0) = (9/5, 12/5); P(4, -3) = 0; P² = P",
     lambda: (lambda P: P == sp.Matrix([[9, 12], [12, 16]]) / 25 and P * vector("5, 0") == vector("9/5, 12/5")
                        and P * vector("4, -3") == sp.zeros(2, 1) and P * P == P)(proyeccion("3,4"))),
    ("reflexión respecto a la recta de (3, 4): H(5, 0) = (-7/5, 24/5) y det H = -1",
     lambda: (lambda H: H * vector("5, 0") == vector("-7/5, 24/5") and H.det() == -1)(reflexion("3,4"))),
    ("núcleo e imagen de [[1, 0, 0], [0, 1, 0]]: rango 2, núcleo (0, 0, 1), dim imagen 2, 1 + 2 = 3",
     lambda: (lambda r: r["rango"] == 2 and len(r["base_nucleo"]) == 1
                        and sp.Matrix.hstack(r["base_nucleo"][0], vector("0, 0, 1")).rank() == 1
                        and r["dim_imagen"] == 2 and r["dim_nucleo"] + r["dim_imagen"] == r["n"] == 3)(
         nucleo_imagen("1,0,0; 0,1,0"))),
    ("T(v) = (v1 + v2, 2) no es lineal: falla en T(0)",
     lambda: (lambda r: not r["lineal"] and "T(0)" in r["falla"])(prueba_linealidad(lambda v: sp.Matrix([v[0] + v[1], 2]), 2))),
    ("T(v) = (v1², v2) no es lineal (falla la aditividad o la homogeneidad)",
     lambda: (lambda r: not r["lineal"] and ("T(u + v)" in r["falla"] or "·u)" in r["falla"]))(
         prueba_linealidad(lambda v: sp.Matrix([v[0] ** 2, v[1]]), 2))),
    ("T(v) = (2v1 - v2, v1 + 3v2): sin contraejemplo, con matriz [[2, -1], [1, 3]]",
     lambda: (lambda r: r["lineal"] and r["matriz"] == sp.Matrix([[2, -1], [1, 3]]))(
         prueba_linealidad(lambda v: sp.Matrix([2 * v[0] - v[1], v[0] + 3 * v[1]]), 2))),
    ("transmisor de 4 a 20 mA, I = 4 + 0.16 T: no es lineal",
     lambda: not prueba_linealidad(lambda v: sp.Matrix([4 + sp.Rational(4, 25) * v[0]]), 1)["lineal"]),
    ("R_z(90 grados)·(1, 0, 0) = (0, 1, 0) y det = 1",
     lambda: (lambda Rz: Rz * vector("1, 0, 0") == vector("0, 1, 0") and Rz.det() == 1)(rotacion3("z", 90, "grados"))),
    ("mano derecha: R_y(90 grados)·(1, 0, 0) = (0, 0, -1) y R_x(90 grados)·(0, 1, 0) = (0, 0, 1)",
     lambda: rotacion3("y", 90, "grados") * vector("1, 0, 0") == vector("0, 0, -1")
             and rotacion3("x", 90, "grados") * vector("0, 1, 0") == vector("0, 0, 1")),
    ("escalamiento(2, 3)·(1, 1) = (2, 3) (k1 en x, k2 en y)",
     lambda: escalamiento(2, 3) * vector("1, 1") == vector("2, 3")),
    ("imagen de [[1, 0, 1], [2, 1, 3], [3, 1, 4]]: columnas ORIGINALES (1, 2, 3) y (0, 1, 1)",
     lambda: nucleo_imagen("1, 0, 1; 2, 1, 3; 3, 1, 4")["base_imagen"] == [vector("1, 2, 3"), vector("0, 1, 1")]),
    ("matriz_desde_imagenes: T(e1) = (2, 1), T(e2) = (-1, 3) → [[2, -1], [1, 3]] (las imágenes son columnas)",
     lambda: matriz_desde_imagenes("2, 1; -1, 3") == sp.Matrix([[2, -1], [1, 3]])),
    ("sin unidad no hay rotación: rotacion(30) y rotacion3(\"z\", 90) → ValueError «indica la unidad»",
     lambda: _mensaje(lambda: rotacion(30)).startswith("indica la unidad")
             and _mensaje(lambda: rotacion3("z", 90)).startswith("indica la unidad")),
    ("ángulos con pi en radianes: «pi/6» = 30 grados, «π/4» = 45 grados, «2*pi/3» = 120 grados; «pi/6» en grados → ValueError",
     lambda: rotacion("pi/6", "radianes") == rotacion(30, "grados") and rotacion("π/4", "radianes") == rotacion(45, "grados")
             and rotacion("2*pi/3", "radianes") == rotacion(120, "grados") and _rechaza(lambda: rotacion("pi/6", "grados"))),
    ("prueba_linealidad con no enteros: parte entera, x + sen(2πx) y x·cos²(πy) no son lineales",
     lambda: all(not prueba_linealidad(f, n)["lineal"] for f, n in _NO_LINEALES)
             and "no se encontró contraejemplo" in prueba_linealidad(lambda v: sp.Matrix([3 * v[0]]), 1)["mensaje"]),
    ("prueba de humo: la calculadora corre con sus valores iniciales, con cada uno de los cinco tipos",
     lambda: all(_humo(calculadora_transformacion, _CONTROLES_4A, tipo=t) for t in _TIPOS)),
    ("prueba de humo: la calculadora de composición corre con sus valores iniciales",
     lambda: _humo(calculadora_composicion, _CONTROLES_4B)),
]
corre_pruebas(PRUEBAS_4, 4)

# %% [markdown]
# **Contrasta (sección 4).** Calcula a mano la imagen de $\mathbf x=(2,0)$ al rotarlo $60^\circ$ en sentido antihorario: $R(60^\circ)\,\mathbf x=(2\cos 60^\circ,\ 2\operatorname{sen} 60^\circ)$. Escribe `mi_imagen = (a, b)` con dos decimales.

# %%
mi_imagen = None     # dos números con dos decimales, por ejemplo (0.25, -1.50)

if mi_imagen is None:     # la referencia se muestra después de tu cálculo
    print("Falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda.")
else:
    try:
        mia = vector(mi_imagen, "tu imagen")
    except ValueError as err:
        print(f"Revisa tu imagen: {str(err).rstrip('.')}.")
    else:
        norma = math.hypot(*[float(e) for e in mia]) if mia.rows == 2 else None
        if mia.rows != 2:          # respuesta imposible: la referencia todavía no se muestra
            print(f"La imagen de un vector del plano tiene dos componentes (escribiste {mia.rows}): revisa tu "
                  "respuesta y vuelve a ejecutar la celda.")
        elif abs(norma - 2) > 0.05:     # una rotación no cambia la longitud: tampoco se muestra la referencia
            print(f"Una rotación no cambia la longitud: la imagen de (2, 0) debe medir 2, y la tuya mide "
                  f"{cifras(norma, 3)}. Revisa tu cálculo y vuelve a ejecutar la celda.")
        else:
            ref = rotacion(60, "grados") * vector("2, 0")
            dos = lambda v: "(" + ", ".join(decimales_fijos(e, 2) for e in v) + ")"
            print(f"La calculadora da: R(60 grados)·(2, 0) = {_tupla(ref)} ≈ {dos(ref)}.")
            con_radianes = rotacion(60, "radianes") * vector("2, 0")
            if all(abs(float(a - b)) <= 0.005 for a, b in zip(mia, ref)):
                print("Coinciden (con tolerancia 0.005).")
            elif all(abs(float(a - b)) <= 0.02 for a, b in zip(mia, con_radianes)):
                print(f"NO coinciden: tu resultado se parece a {dos(con_radianes)}, que sale si la calculadora "
                      "está en radianes (cos 60 rad ≈ -0.952). Probablemente usaste radianes: pon tu calculadora en grados.")
            else:
                print("NO coinciden: cos 60° = 1/2 y sen 60° = √3/2 ≈ 0.866; multiplica cada uno por 2.")
            print("\nEl procedimiento completo:")
            _informe_transformacion(rotacion(60, "grados"), "R(60 grados)", vector("2, 0"), 4, tex=r"R(60^\circ)")
