# ID: ALG-U6-NB06
# Notebook: alg/u6_aplicaciones.ipynb · sección 6.6 transformaciones de coordenadas
# Repositorio: alg/u6_aplicaciones/06_transformaciones.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6. Transformaciones de coordenadas
#
# **Un marco de referencia** (el de una cámara, el de la pinza de un robot) se describe respecto a otro, fijo, con dos cosas:
#
# - una **rotación** $R$ de $3\times3$, cuyas columnas son los ejes $x$, $y$, $z$ del marco móvil escritos en el fijo. Por eso son ortonormales: $R^TR=I$, y además $\det R=1$ (si fuera $-1$, sería una reflexión);
# - una **traslación** $\mathbf t$: dónde está el origen del marco móvil, escrito en el fijo.
#
# Un punto que en el marco móvil $C$ tiene coordenadas $\mathbf p_C$ tiene en el fijo $B$ las coordenadas
#
# $$\mathbf p_B=R\,\mathbf p_C+\mathbf t .$$
#
# Las rotaciones básicas son antihorarias vistas desde la punta de su eje (regla de la mano derecha); por ejemplo, $R_z(\theta)=\begin{bmatrix}\cos\theta&-\operatorname{sen}\theta&0\\ \operatorname{sen}\theta&\cos\theta&0\\0&0&1\end{bmatrix}$.
#
# **Coordenadas homogéneas.** Para que rotar y trasladar sea **un solo producto de matrices**, se usa la matriz de $4\times4$
#
# $$T=\begin{bmatrix}R&\mathbf t\\ \mathbf 0^T&1\end{bmatrix},\qquad T\begin{bmatrix}\mathbf p\\1\end{bmatrix}=\begin{bmatrix}R\,\mathbf p+\mathbf t\\1\end{bmatrix},\qquad T\begin{bmatrix}\mathbf d\\0\end{bmatrix}=\begin{bmatrix}R\,\mathbf d\\0\end{bmatrix}.$$
#
# Un **punto** lleva un 1 al final y se traslada; una **dirección** (una velocidad, el eje de una herramienta) lleva un 0 y solo se rota.
#
# **Composición.** Si $T_{AB}$ describe al marco $B$ en el $A$ y $T_{BC}$ al $C$ en el $B$, entonces $T_{AC}=T_{AB}\,T_{BC}$. Aplicada a un punto, actúa primero la de la derecha. **El orden importa:** en general $T_1T_2\neq T_2T_1$.
#
# **Inversa.** Para regresar del marco fijo al móvil:
#
# $$T^{-1}=\begin{bmatrix}R^T&-R^T\mathbf t\\ \mathbf 0^T&1\end{bmatrix}.$$
#
# Ojo: **no es la transpuesta** de $T$. Solo la parte de rotación se invierte transponiendo; la traslación cambia a $-R^T\mathbf t$. La comprobación es $T\,T^{-1}=I$, que conviene hacer con una rotación y una traslación que no sean triviales.

# %%
_UNIDADES = {"grados": "grados", "grado": "grados", "radianes": "radianes", "radian": "radianes", "rad": "radianes"}
_AYUDA_ANGULO = "escribe un número (30, 37.5, 1/3) o, en radianes, un múltiplo de pi (pi/6, 2*pi/3, π/4, -pi)"


def _unidad(unidad):
    """'grados' o 'radianes' (acepta mayúsculas, acentos, el singular y 'rad'). No hay unidad por omisión: sin unidad,
    o con otra cosa, ValueError en español."""
    if unidad is None or (isinstance(unidad, str) and not unidad.strip()):
        raise ValueError("indica la unidad: 'grados' o 'radianes' (no hay unidad por omisión)")
    t = normaliza_texto(unidad) if isinstance(unidad, str) else ""
    if t not in _UNIDADES:
        raise ValueError(f"indica la unidad: 'grados' o 'radianes' (escribiste «{unidad}»)")
    return _UNIDADES[t]


def _limpia_angulo(texto):
    """Quita el signo de grados (90°), los paréntesis que envuelven todo ((pi/2)) y lee '·' como producto (2·pi/3)."""
    t = str(texto).replace("°", "").replace("º", "").replace("·", "*").replace("×", "*").strip()
    while t.startswith("(") and t.endswith(")"):
        dentro, nivel = t[1:-1], 0
        for ch in dentro:                       # solo si los paréntesis de afuera son pareja: '(1)/(2)' no se toca
            nivel += 1 if ch == "(" else (-1 if ch == ")" else 0)
            if nivel < 0:
                break
        if nivel != 0:
            break
        t = dentro.strip()
    return t


def angulo_exacto(angulo):
    """El ángulo como número exacto: 30 -> 30, '37.5' -> 75/2, '90°' -> 90, 'pi/6' -> π/6, '2*pi/3', '2·pi/3', '2π/3',
    '(pi/2)', '-π/4', sp.pi/2. El lector es restringido: solo números, pi, productos y cocientes (no evalúa expresiones
    libres)."""
    if isinstance(angulo, sp.Basic) and not angulo.is_Rational and not angulo.is_Float:
        if angulo.is_number and (angulo / sp.pi).is_Rational:
            return angulo
        raise ValueError(f"«{angulo}» no es un ángulo que entienda esta calculadora: {_AYUDA_ANGULO}")
    if not isinstance(angulo, str):
        return numero(angulo, "el ángulo")
    limpio = _limpia_angulo(angulo)
    if "(" in limpio or ")" in limpio:
        raise ValueError(f"«{angulo}» no es un ángulo que entienda esta calculadora: {_AYUDA_ANGULO}")
    t = _normaliza_entrada(limpio).replace(" ", "").lower().replace("π", "pi")
    if "pi" not in t:
        try:
            return numero(limpio, "el ángulo")
        except ValueError:
            raise ValueError(f"«{angulo.strip()}» no es un ángulo que entienda esta calculadora: {_AYUDA_ANGULO}") from None
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
    if u == "radianes" and isinstance(angulo, str) and ("°" in angulo or "º" in angulo):
        raise ValueError("escribiste el ángulo con °, pero la unidad es 'radianes': cambia la unidad a 'grados' o quita "
                         "el °")
    a = angulo_exacto(angulo)
    if u == "grados":
        if a.has(sp.pi):
            raise ValueError("con la unidad 'grados' escribe el ángulo sin pi (por ejemplo 30); pi/6 ya está en radianes")
        return sp.rad(a)
    return a


def _texto_angulo(angulo, unidad):
    """'90 grados', 'π/2 radianes'."""
    a = angulo_exacto(angulo)
    return f"{_ntxt(a) if a.is_Rational else _txt(a).replace('·π', 'π')} {_unidad(unidad)}"


def rot(eje, angulo, unidad=None):
    """Rotación exacta de R³ alrededor del eje 'x', 'y' o 'z', antihoraria vista desde la punta del eje (regla de la
    mano derecha). La unidad es OBLIGATORIA: rot('z', 90, 'grados') o rot('z', 'pi/2', 'radianes'); sin unidad,
    ValueError «indica la unidad». Con ángulos notables (30, 45, 60, 90 grados...) las entradas son exactas."""
    u = _unidad(unidad)
    e = normaliza_texto(eje) if isinstance(eje, str) else ""
    if e not in ("x", "y", "z"):
        raise ValueError(f"el eje de rotación debe ser 'x', 'y' o 'z' (escribiste «{eje}»)")
    t = _angulo_rad(angulo, u)
    c, s = sp.cos(t), sp.sin(t)
    if sp.count_ops(c) + sp.count_ops(s) > 20:     # sympy los desarrolló en radicales anidados (1.5°, 12°...): se dejan
        c, s = sp.cos(t, evaluate=False), sp.sin(t, evaluate=False)          # como cos(π/120), que se calcula rápido
    if e == "x":
        return sp.Matrix([[1, 0, 0], [0, c, -s], [0, s, c]])
    if e == "y":
        return sp.Matrix([[c, 0, s], [0, 1, 0], [-s, 0, c]])
    return sp.Matrix([[c, -s, 0], [s, c, 0], [0, 0, 1]])


def _matriz_exacta(M, nombre):
    """Una Matrix de sympy se toma tal cual (puede traer senos, cosenos y raíces exactos); texto o listas se leen como
    racionales exactos."""
    if isinstance(M, sp.MatrixBase):
        if not all(sp.sympify(e).is_number for e in M):
            raise ValueError(f"{nombre} debe tener solo números")
        return sp.Matrix(M)
    return como_matriz(M, nombre)


def _motivo_no_rotacion(R):
    """None si R es una rotación (RᵀR = I y det R = 1); si no, el motivo."""
    if R.shape != (3, 3):
        return f"debe ser de 3×3 (es de {dims(R)})"
    if not all(_es_cero(e) for e in R.T * R - sp.eye(3)):
        return ("RᵀR ≠ I: sus columnas no son ortonormales, así que no conserva longitudes ni ángulos (si escribiste "
                "decimales redondeados, construye R con rot())")
    d = R.det()
    if _es_cero(d + 1):
        return "det R = -1: es una reflexión (cambia la orientación), no una rotación"
    return None if _es_cero(d - 1) else f"det R = {_valor(d)}, no 1"


def homogenea(R, t):
    """La matriz homogénea T = [[R, t], [0, 0, 0, 1]] de 4×4. R debe ser una rotación (RᵀR = I y det R = 1); si no,
    ValueError «R no es una rotación: ...». t es la traslación, un vector de R³."""
    R = _matriz_exacta(R, "R")
    motivo = _motivo_no_rotacion(R)
    if motivo:
        raise ValueError(f"R no es una rotación: {motivo}")
    t = _vector3(t, "la traslación t")
    T = sp.eye(4)
    T[:3, :3] = R
    T[:3, 3] = t
    return T


def _valida_T(T, nombre="T"):
    """T como Matrix de 4×4 homogénea: último renglón (0, 0, 0, 1) y bloque R de rotación."""
    T = _matriz_exacta(T, nombre)
    if T.shape != (4, 4):
        raise ValueError(f"{nombre} debe ser de 4×4 (es de {dims(T)})")
    if list(T.row(3)) != [0, 0, 0, 1]:
        raise ValueError(f"{nombre} no es homogénea: su último renglón debe ser (0, 0, 0, 1)")
    motivo = _motivo_no_rotacion(T[:3, :3])
    if motivo:
        raise ValueError(f"el bloque R de {nombre} no es una rotación: {motivo}")
    return T


def _tipo(tipo):
    t = normaliza_texto(tipo) if isinstance(tipo, str) else ""
    if t not in ("punto", "direccion"):
        raise ValueError("indica el tipo: 'punto' (lleva 1 al final y se traslada) o 'dirección' (lleva 0 y solo se rota)")
    return t


def transforma(T, p, tipo=None):
    """Aplica T a p. El tipo es OBLIGATORIO: 'punto' (p se escribe [p; 1] y se traslada) o 'dirección' ([p; 0], solo
    se rota). Devuelve el vector de R³ resultante, exacto."""
    T, tp = _valida_T(T), _tipo(tipo)
    p = _vector3(p, "el punto" if tp == "punto" else "la dirección")
    h = p.col_join(sp.Matrix([1 if tp == "punto" else 0]))
    return _simplifica(T * h)[:3, 0]


def inversa_homogenea(T):
    """Inversa de T con la fórmula T⁻¹ = [[Rᵀ, -Rᵀt], [0, 0, 0, 1]], sin invertir una matriz de 4×4. Diccionario con
    Tinv, R, t, traslacion (-Rᵀt), comprobacion (T·T⁻¹, simplificada) y es_identidad (True si T·T⁻¹ = I)."""
    T = _valida_T(T)
    R, t = T[:3, :3], T[:3, 3]
    Ti = sp.eye(4)
    Ti[:3, :3] = R.T
    Ti[:3, 3] = _simplifica(-R.T * t)
    C = _simplifica(T * Ti)
    return {"Tinv": Ti, "R": R, "t": t, "traslacion": Ti[:3, 3], "comprobacion": C,
            "es_identidad": all(_es_cero(x) for x in C - sp.eye(4))}


def compone(*Ts, nombres=None):
    """Producto T1·T2·...·Tk en el orden dado. Diccionario con T (el producto, exacto) y texto (en qué orden actúan
    sobre un punto: primero la de la derecha)."""
    if not Ts:
        raise ValueError("escribe al menos una transformación")
    if len(Ts) > 6:
        raise ValueError(f"escribiste {len(Ts)} transformaciones; esta función compone hasta 6")
    Ms = [_valida_T(T, f"T{k}") for k, T in enumerate(Ts, 1)]
    nombres = list(nombres) if nombres is not None else [f"T{k}" for k in range(1, len(Ms) + 1)]
    if len(nombres) != len(Ms):
        raise ValueError(f"hay {_cuenta(len(Ms), 'transformación', 'transformaciones')} y "
                         f"{_cuenta(len(nombres), 'nombre', 'nombres')}: debe haber un nombre por cada una")
    P = Ms[0]
    for M in Ms[1:]:
        P = _simplifica(P * M)
    if len(Ms) == 1:
        texto = f"{nombres[0]}: una sola transformación"
    else:
        orden = list(reversed(nombres))
        if len(orden) == 2:
            texto = f"{'·'.join(nombres)}: aplicada a un punto, actúa primero {orden[0]} y después {orden[1]} (el orden importa)"
        else:
            texto = (f"{'·'.join(nombres)}: aplicada a un punto, actúa primero {orden[0]}"
                     + "".join(f", después {nm}" for nm in orden[1:-1]) + f" y al final {orden[-1]} (el orden importa)")
    return {"T": P, "texto": texto}


def _informe_transformacion(R, T, p, tipo, angulo_texto, eje, n_cifras):
    tp = _tipo(tipo)
    t = T[:3, 3]
    muestra(f"R = R_{eje}({angulo_texto})", R, n_cifras, tex=rf"R=R_{{{eje}}}(\text{{{angulo_texto}}})")
    print(f"t = {_tupla_u(t, 'm', n_cifras)}  (el origen del marco móvil, en el marco fijo)")
    muestra("T", T, n_cifras, tex=r"T=\begin{bmatrix}R&\mathbf t\\ \mathbf 0^T&1\end{bmatrix}")
    q = transforma(T, p, tp)
    if tp == "punto":
        print(f"Punto p = {_tupla_u(p, 'm')} (en el marco móvil) → [p; 1] → T·[p; 1] da R·p + t = "
              f"{_tupla_u(q, 'm', n_cifras)} (en el marco fijo)")
    else:
        print(f"Dirección d = {_tupla(p)} → [d; 0] → T·[d; 0] da R·d = {_tupla(q, n_cifras)}: la traslación no la "
              "mueve, solo se rota")
    inv = inversa_homogenea(T)
    muestra("T⁻¹", inv["Tinv"], n_cifras, tex=r"T^{-1}=\begin{bmatrix}R^T&-R^T\mathbf t\\ \mathbf 0^T&1\end{bmatrix}")
    print(f"-Rᵀt = {_tupla_u(inv['traslacion'], 'm', n_cifras)};  comprobación: T·T⁻¹ " +
          ("= I" if inv["es_identidad"] else "≠ I (¡revisa!)"))
    if not _es_nulo(t):
        print("Tᵀ no es la inversa: T·Tᵀ ≠ I, porque t ≠ 0.")
    else:
        print("Con t = 0 solo hay rotación, y entonces Tᵀ sí es la inversa: T·Tᵀ = I.")
    return q


def _redondo(o):
    """Un origen en decimales como Matrix con racionales cortos (2.0 -> 2), solo para rotularlo."""
    return sp.Matrix([sp.Rational(repr(float(v))) if abs(v) < 1e15 else v for v in o])


def dibuja_marcos(marcos, puntos=None, n_cifras=4, unidad="m"):
    """Dibujo 3D de hasta 6 marcos (nombre, T): sus ejes x, y, z (continuo, discontinuo y punteado), con un color y
    un marcador por marco. puntos: lista opcional de (rótulo, vector, 'punto' o 'dirección') en el marco fijo; un
    punto se marca con una estrella y una dirección con una flecha desde el origen."""
    if len(marcos) > 6:
        print("No se dibuja: son más de 6 marcos.")
        return
    datos_m = []
    for nombre, T in marcos:
        T = _valida_T(T)
        datos_m.append((nombre, _flotante(T[:3, :3]), _flotante(T[:3, 3]).ravel()))
    origenes = np.array([o for _, _, o in datos_m])
    extension = float(np.ptp(origenes, axis=0).max()) if len(datos_m) > 1 else 0.0
    largo = max(1.0, 0.45 * extension)
    fig = plt.figure(figsize=(7.4, 8.2))
    ax = fig.add_subplot(projection="3d")
    fig.subplots_adjust(bottom=0.2)
    todos = [origenes]
    estilos_eje = ("-", "--", ":")
    for k, (nombre, R, o) in enumerate(datos_m):
        color = _COLORES[k % len(_COLORES)]
        for i, letra in enumerate("xyz"):
            fin = o + (0.8 if k % 2 == 0 else 0.6) * largo * R[:, i]   # largos distintos: las puntas no se enciman
            ax.plot(*np.array([o, fin]).T, linestyle=estilos_eje[i], color=color, lw=2.4)
            ax.text(*fin, f" {letra}{_sub(k)}", fontsize=9, color="black")
            todos.append(fin[None, :])
        ax.scatter(*o, marker=_MARCAS[k % len(_MARCAS)], s=60, color=color, depthshade=False,
                   label=f"marco {nombre}: origen {_tupla_dec(_redondo(o), n_cifras)} {unidad}")
    for etiqueta, v, tp in (puntos or []):
        vf = _flotante(sp.Matrix(v)).ravel()
        if normaliza_texto(tp) == "punto":
            ax.scatter(*vf, marker="*", s=120, color="black", depthshade=False, zorder=5,
                       label=f"{etiqueta}: {_tupla_dec(sp.Matrix(v), n_cifras)} {unidad}")
            todos.append(vf[None, :])
        else:
            punta = largo * vf / max(np.linalg.norm(vf), 1e-300)
            ax.quiver(0, 0, 0, *punta, color="black", linewidth=2.6, arrow_length_ratio=0.22)
            ax.plot([], [], "-", color="black", lw=2.6, label=f"{etiqueta}: {_tupla_dec(sp.Matrix(v), n_cifras)} (flecha)")
            todos.append(punta[None, :])
    for ls, letra in zip(estilos_eje, "xyz"):
        ax.plot([], [], linestyle=ls, color="0.3", lw=2.4, label=f"eje {letra} de cada marco")
    P = np.vstack(todos)
    centro = (P.max(axis=0) + P.min(axis=0)) / 2
    radio = max(float((P.max(axis=0) - P.min(axis=0)).max()) / 2, 1.0) * 1.15
    ax.set_xlim(centro[0] - radio, centro[0] + radio)
    ax.set_ylim(centro[1] - radio, centro[1] + radio)
    ax.set_zlim(centro[2] - radio, centro[2] + radio)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xlabel(f"x ({unidad})"); ax.set_ylabel(f"y ({unidad})"); ax.set_zlabel(f"z ({unidad})")
    ax.view_init(elev=25, azim=-58)
    ax.set_title("Marcos de referencia (ejes escritos en el marco fijo)", fontsize=10)
    fig.legend(*ax.get_legend_handles_labels(), loc="lower center", ncol=2, fontsize=8)
    plt.show()

# %% [markdown]
# **Recorrido paso a paso.** (1) Un marco móvil girado $90^\circ$ alrededor de $z$ y con el origen en $\mathbf t=(2,0,0)$: el punto $(1,0,0)$ del marco móvil está en $(2,1,0)$ en el fijo (se rota a $(0,1,0)$ y se traslada), y la dirección $(1,0,0)$ pasa a $(0,1,0)$: se rota, pero no se traslada. (2) La inversa de $T$ con $R=R_z(90^\circ)$ y $\mathbf t=(1,2,3)$: su traslación es $-R^T\mathbf t=(-2,\,1,\,-3)$, y $T\,T^{-1}=I$; en cambio, $T\,T^T\neq I$. (3) Dos composiciones con $T_1$ = rotación $R_x(90^\circ)$ y $T_2$ = traslación $(0,0,1)$, aplicadas al punto $(0,1,0)$: $T_1T_2$ primero traslada y después rota, y $T_2T_1$ primero rota y después traslada. Dan puntos distintos.

# %%
print("(1) R = R_z(90 grados), t = (2, 0, 0)")
_R90 = rot("z", 90, "grados")
_T1 = homogenea(_R90, (2, 0, 0))
_q1 = _informe_transformacion(_R90, _T1, vector("1, 0, 0"), "punto", "90 grados", "z", 4)
print(f"La dirección (1, 0, 0), en cambio: {_tupla(transforma(_T1, (1, 0, 0), 'dirección'))}")
dibuja_marcos([("0 (fijo)", sp.eye(4)), ("1 (móvil)", _T1)], puntos=[("p en el marco fijo", _q1, "punto")])
print("=" * 90 + "\n\n(2) Inversa con R = R_z(90 grados) y t = (1, 2, 3)")
_T2 = homogenea(_R90, (1, 2, 3))
_inv2 = inversa_homogenea(_T2)
muestra("T", _T2)
muestra("T⁻¹", _inv2["Tinv"], tex=r"T^{-1}")
print(f"-Rᵀt = {_tupla(_inv2['traslacion'])};  T·T⁻¹ = I: {'sí' if _inv2['es_identidad'] else 'no'};  "
      f"T·Tᵀ = I: {'sí' if _T2 * _T2.T == sp.eye(4) else 'no'}")
muestra("T·Tᵀ", _T2 * _T2.T, tex=r"T\,T^T")
print("=" * 90 + "\n\n(3) El orden importa: T1 = rotación R_x(90 grados), T2 = traslación (0, 0, 1), punto (0, 1, 0)")
_Rx = homogenea(rot("x", 90, "grados"), (0, 0, 0))
_Tz = homogenea(sp.eye(3), (0, 0, 1))
for _orden in ((_Rx, _Tz), (_Tz, _Rx)):
    _c = compone(*_orden, nombres=["T1", "T2"] if _orden[0] is _Rx else ["T2", "T1"])
    print(f"{_c['texto']}  →  {_tupla(transforma(_c['T'], (0, 1, 0), 'punto'))}")

# %% [markdown]
# **Calculadora.** Elige el eje y el ángulo de la rotación **con su unidad** (en radianes puedes escribir `pi/2` o `π/4`), la traslación $\mathbf t$ en metros, el vector y si es un **punto** o una **dirección**. La calculadora arma $T$, transforma el vector, calcula $T^{-1}$ con la fórmula y comprueba $T\,T^{-1}=I$, y dibuja los dos marcos. Prueba 90 en radianes (no es lo mismo que 90 grados), cambiar el punto por una dirección y una traslación $(0,0,0)$, con la que $T^T$ sí es la inversa.

# %%
def calculadora_transformacion(eje, angulo_txt, unidad, traslacion_txt, punto_txt, tipo, n_cifras):
    with bloque():                 # el texto aparece completo, de una sola vez; el dibujo va después
        try:
            R = rot(eje, angulo_txt, unidad)
            T = homogenea(R, traslacion_txt)
            p = _vector3(punto_txt, "el punto" if _tipo(tipo) == "punto" else "la dirección")
            texto = _texto_angulo(angulo_txt, unidad)
        except ValueError as err:
            _reporta(err); return
        _aviso_coma(traslacion_txt, punto_txt)
        q = _informe_transformacion(R, T, p, tipo, texto, normaliza_texto(eje), n_cifras)
    dibuja_marcos([("0 (fijo)", sp.eye(4)), ("1 (móvil)", T)],
                  puntos=[("p en el marco fijo" if _tipo(tipo) == "punto" else "R·d", q, tipo)], n_cifras=n_cifras)


_CONTROLES_6 = dict(
    eje=widgets.Dropdown(options=["x", "y", "z"], value="z", description="eje:"),
    angulo_txt=widgets.Text(value="90", description="ángulo:", continuous_update=False),
    unidad=widgets.ToggleButtons(options=["grados", "radianes"], value="grados", description="unidad:"),
    traslacion_txt=widgets.Text(value="2, 0, 0", description="t (m) =", continuous_update=False),
    punto_txt=widgets.Text(value="1, 0, 0", description="punto/dirección:", continuous_update=False,
                           style={"description_width": "initial"}),
    tipo=widgets.ToggleButtons(options=["punto", "dirección"], value="punto", description="es un:"),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_transformacion, **_CONTROLES_6);

# %%
# Casos de prueba de la sección 6 (resultado conocido)
def _rapida(angulo, tope=6.0):
    """La calculadora con un ángulo cuyo coseno sympy desarrolla en radicales anidados: termina en menos de tope
    segundos, sin «Revisa la entrada», y comprueba T·T⁻¹ = I."""
    inicio = time.time()
    texto, error = _captura(calculadora_transformacion, "z", angulo, "grados", "1, 2, 3", "1, 1, 1", "punto", 4)
    return (error is None and time.time() - inicio < tope and "Revisa la entrada" not in texto
            and "comprobación: T·T⁻¹ = I" in texto)


PRUEBAS_6 = [
    ("rot('z', 90, 'grados') = [[0, -1, 0], [1, 0, 0], [0, 0, 1]], y en radianes con pi/2 (y con 'pi/2', '(pi/2)' o "
     "'90°' en grados) es igual; R_y(90 grados)·(0, 0, 1) = (1, 0, 0); '2·pi/3' = 2π/3",
     lambda: rot("z", 90, "grados") == sp.Matrix([[0, -1, 0], [1, 0, 0], [0, 0, 1]])
             and rot("z", sp.pi / 2, "radianes") == rot("z", 90, "grados") == rot("z", "pi/2", "radianes")
             and rot("z", "(pi/2)", "radianes") == rot("z", "90°", "grados") == rot("z", 90, "grados")
             and rot("y", 90, "grados") * vector("0, 0, 1") == vector("1, 0, 0")
             and angulo_exacto("2·pi/3") == 2 * sp.pi / 3),
    ("R_z(90 grados) con t = (2, 0, 0): el punto (1, 0, 0) → (2, 1, 0) y la dirección (1, 0, 0) → (0, 1, 0)",
     lambda: transforma(_T1, (1, 0, 0), "punto") == vector("2, 1, 0")
             and transforma(_T1, (1, 0, 0), "dirección") == vector("0, 1, 0")),
    ("inversa_homogenea con R_z(90 grados) y t = (1, 2, 3) → traslación (-2, 1, -3) y T·T⁻¹ = I",
     lambda: _inv2["traslacion"] == vector("-2, 1, -3") and _inv2["es_identidad"]
             and _inv2["Tinv"] == _T2.inv()),
    ("la transpuesta de T no es su inversa cuando t ≠ 0: T·Tᵀ ≠ I (y con t = 0 sí lo es)",
     lambda: _T2 * _T2.T != sp.eye(4) and (lambda T: T * T.T == sp.eye(4))(homogenea(_R90, (0, 0, 0)))),
    ("composición: T1·T2 ≠ T2·T1 con una rotación y una traslación no nulas, y (T1·T2)⁻¹ = T2⁻¹·T1⁻¹",
     lambda: compone(_Rx, _Tz)["T"] != compone(_Tz, _Rx)["T"]
             and inversa_homogenea(compone(_Rx, _Tz)["T"])["Tinv"] ==
                 _simplifica(inversa_homogenea(_Tz)["Tinv"] * inversa_homogenea(_Rx)["Tinv"])
             and transforma(compone(_Rx, _Tz)["T"], (0, 1, 0), "punto") == vector("0, -1, 1")
             and transforma(compone(_Tz, _Rx)["T"], (0, 1, 0), "punto") == vector("0, 0, 2")),
    ("homogenea con R = [[1, 1, 0], [0, 1, 0], [0, 0, 1]] → ValueError; con la reflexión diag(1, 1, -1) → ValueError "
     "(det = -1)",
     lambda: _mensaje(lambda: homogenea([[1, 1, 0], [0, 1, 0], [0, 0, 1]], (0, 0, 0))).startswith("R no es una rotación")
             and (lambda m: m.startswith("R no es una rotación") and "det R = -1" in m and "reflexión" in m)(
                 _mensaje(lambda: homogenea(sp.diag(1, 1, -1), (0, 0, 0))))),
    ("rot('z', 90) sin unidad → ValueError «indica la unidad»; con 'sexagesimales', también",
     lambda: _mensaje(lambda: rot("z", 90)).startswith("indica la unidad")
             and _mensaje(lambda: rot("z", 90, "sexagesimales")).startswith("indica la unidad")),
    ("prueba de humo: la calculadora corre con sus valores iniciales, con una dirección, en radianes y con 1.5 y 12 "
     "grados (cada uno en menos de 6 s y con T·T⁻¹ = I), sin error y sin «Revisa la entrada»",
     lambda: _humo(calculadora_transformacion, _CONTROLES_6)
             and _humo(calculadora_transformacion, _CONTROLES_6, tipo="dirección", eje="x", n_cifras=6)
             and _humo(calculadora_transformacion, _CONTROLES_6, unidad="radianes", angulo_txt="pi/3")
             and all(_rapida(angulo) for angulo in ("1.5", "12"))),
]
corre_pruebas(PRUEBAS_6, 6)

# %% [markdown]
# **Contrasta (sección 6).** Con $R=R_z(90^\circ)$ y $\mathbf t=(0,0,2)$, transforma a mano el punto $(1,1,0)$: primero $R\,\mathbf p$, después suma $\mathbf t$. Escribe `mi_punto = (x, y, z)`.

# %%
mi_punto = None     # tres números (x, y, z), por ejemplo (0.5, -1, 3)

if mi_punto is None:     # la referencia se muestra después de tu cálculo
    print("Falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda.")
else:
    try:
        mio = vector(mi_punto, "mi_punto")
    except ValueError as err:
        print(f"Revisa tu respuesta: {str(err).rstrip('.')}.")
    else:
        if mio.rows != 3:          # respuesta imposible: la referencia todavía no se muestra
            print(f"El punto transformado es de R³, con tres coordenadas (escribiste {mio.rows}): revisa tu respuesta y "
                  "vuelve a ejecutar la celda.")
        else:
            T_c = homogenea(rot("z", 90, "grados"), (0, 0, 2))          # los datos del enunciado
            ref = transforma(T_c, (1, 1, 0), "punto")
            print(f"La calculadora da: R·p + t = {_tupla(ref)}.")
            al_reves = transforma(homogenea(rot("z", -90, "grados"), (0, 0, 2)), (1, 1, 0), "punto")
            if all(abs(float(a - b)) <= 0.01 for a, b in zip(mio, ref)):
                print("Coinciden (con tolerancia de 0.01).")
            elif all(abs(float(a - b)) <= 0.01 for a, b in zip(mio, al_reves)):
                print("NO coinciden: giraste en sentido horario. Revisa el signo de la rotación: R_z(90°) gira en sentido "
                      "antihorario visto desde +z, así que manda (1, 0, 0) a (0, 1, 0) y (0, 1, 0) a (-1, 0, 0).")
            elif all(abs(float(a - b)) <= 0.01 for a, b in zip(mio, ref - vector("0, 0, 2"))):
                print("NO coinciden: te faltó sumar la traslación t = (0, 0, 2).")
            else:
                print("NO coinciden: calcula R·p con R = [[0, -1, 0], [1, 0, 0], [0, 0, 1]] y después suma t.")
            print("\nEl procedimiento completo:")
            _informe_transformacion(rot("z", 90, "grados"), T_c, vector("1, 1, 0"), "punto", "90 grados", "z", 4)
