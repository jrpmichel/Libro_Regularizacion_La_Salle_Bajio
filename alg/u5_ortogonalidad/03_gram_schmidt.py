# ID: ALG-U5-NB03
# Notebook: alg/u5_ortogonalidad.ipynb · sección 5.3 proceso de Gram-Schmidt
# Repositorio: alg/u5_ortogonalidad/03_gram_schmidt.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Proceso de Gram-Schmidt
#
# **Idea.** A partir de vectores $\mathbf v_1,\dots,\mathbf v_k$ se construyen vectores **ortogonales** $\mathbf u_1,\dots,\mathbf u_k$ que generan lo mismo: a cada $\mathbf v_j$ se le quita su proyección sobre los $\mathbf u$ anteriores (que ya son ortogonales, así que basta sumar proyecciones, sección 2) y queda solo su parte perpendicular:
#
# $$\mathbf u_1=\mathbf v_1,\qquad \mathbf u_j=\mathbf v_j-\sum_{i<j}\frac{\mathbf v_j\cdot\mathbf u_i}{\mathbf u_i\cdot\mathbf u_i}\,\mathbf u_i,\qquad \mathbf q_j=\frac{\mathbf u_j}{\|\mathbf u_j\|}.$$
#
# - En cada paso, $\operatorname{gen}\{\mathbf u_1,\dots,\mathbf u_j\}=\operatorname{gen}\{\mathbf v_1,\dots,\mathbf v_j\}$, porque $\mathbf u_j$ es $\mathbf v_j$ menos una combinación de los anteriores.
# - Si $\mathbf v_j$ es combinación de los anteriores, su parte perpendicular es cero: $\mathbf u_j=\mathbf 0$, se descarta y se sigue con el siguiente. Así, el proceso también detecta dependencia.
# - El resultado depende del **orden** de los vectores (el primero se queda como está), pero el subespacio que generan no.
# - Los $\mathbf u_j$ son ortogonales; los $\mathbf q_j$ son, además, de longitud 1 (ortonormales). Para no arrastrar raíces en las cuentas, conviene normalizar al final.
#
# **Factorización QR.** Si las columnas de $A=[\,\mathbf v_1\ \cdots\ \mathbf v_k\,]$ son independientes, $A=QR$ con
#
# $$Q=[\,\mathbf q_1\ \cdots\ \mathbf q_k\,]\quad(Q^TQ=I),\qquad R=Q^TA,\quad r_{ij}=\mathbf q_i\cdot\mathbf v_j .$$
#
# $R$ es triangular superior: $r_{ij}=0$ si $i>j$, porque $\mathbf v_j$ está en $\operatorname{gen}\{\mathbf q_1,\dots,\mathbf q_j\}$, que es perpendicular a $\mathbf q_i$. Su diagonal es positiva, $r_{jj}=\|\mathbf u_j\|$, así que $R$ es invertible. Si las columnas son dependientes, algún $\mathbf u_j$ es cero y no hay QR con $R$ invertible.
#
# **Nota numérica.** En decimales, Gram-Schmidt clásico pierde ortogonalidad con vectores casi dependientes: al restar dos vectores casi iguales, el error de redondeo queda del tamaño del resultado. Por eso numpy no lo usa: `np.linalg.qr` trabaja con reflexiones de Householder, que conservan la ortogonalidad, y puede devolver columnas de $Q$ y renglones de $R$ con el signo cambiado (la diagonal de su $R$ puede ser negativa). Cambiar el signo de una columna de $Q$ y del renglón correspondiente de $R$ no cambia el producto $QR$.

# %%
def _resta(c, nombre):
    """El término que se resta, con el signo que toca: c = 1/2 -> '- (1/2)·u1'; c = -2 -> '+ 2·u1'; c = 0 -> '- 0·u1'."""
    return f"- 0·{nombre}" if c == 0 else _termino(-c, nombre, False)


def gram_schmidt(vectores, normalizar=False, n_cifras=4):
    """Gram-Schmidt exacto. Diccionario con u (los u_j distintos de cero, exactos), indices (el j de cada uno: u_j
    viene de v_j), q (los u_j normalizados, exactos, con radicales), normas (‖u_j‖), pasos (textos legibles:
    'u2 = v2 - (1/2)·u1 = (1/2, -1/2, 1)'; con normalizar=True, también 'q1 = u1/‖u1‖ = ...'), dependientes (los j
    de los v_j que dieron u_j = 0, desde 1), coeficientes ({j: [(i, (v_j·u_i)/(u_i·u_i)), ...]}), gram (matriz de
    Gram de los u, diagonal) y A."""
    A = lista_vectores(vectores)
    n, k = A.shape
    us, indices, pasos, dependientes, coefs = [], [], [], [], {}
    for j in range(1, k + 1):
        v = A[:, j - 1]
        cs = [_punto(v, u) / _punto(u, u) for u in us]
        w = v - sum((c * u for c, u in zip(cs, us)), sp.zeros(n, 1))
        coefs[j] = list(zip(indices, cs))
        texto = f"u{j} = v{j}" + "".join(" " + _resta(c, f"u{i}") for i, c in zip(indices, cs)) + \
                f" = {_tupla(w, None if _exacta_legible(w) else n_cifras)}"
        if _es_nulo(w):
            dependientes.append(j)
            texto += (f": v{j} es el vector cero y se descarta" if _es_nulo(v) else
                      f": v{j} es combinación de los anteriores y se descarta")
        else:
            us.append(w)
            indices.append(j)
        pasos.append(texto)
    normas = [sp.sqrt(_punto(u, u)) for u in us]
    q = [u / nu for u, nu in zip(us, normas)]
    if normalizar:
        for i, u, nu, qq in zip(indices, us, normas, q):
            pasos.append(f"q{i} = u{i}/‖u{i}‖, con ‖u{i}‖ = {_valor(nu, n_cifras)}:  q{i} = {_tupla(qq, n_cifras)}")
    U = sp.Matrix.hstack(*us) if us else sp.zeros(n, 0)
    return {"u": us, "indices": indices, "q": q, "normas": normas, "pasos": pasos, "dependientes": dependientes,
            "coeficientes": coefs, "gram": U.T * U, "A": A}


def columnas_ortonormales(Q):
    """True si las columnas de Q son ortonormales, es decir, si QᵀQ = I (comprobado en forma exacta)."""
    Q = sp.Matrix(Q)
    return all(_es_cero(e) for e in Q.T * Q - sp.eye(Q.cols))


def factorizacion_qr(vectores):
    """A = QR exacta, con los vectores como columnas de A: Q = [q1 ... qk] (columnas ortonormales) y R = QᵀA
    (triangular superior, diagonal positiva). Diccionario con Q, R, A, QtQ_es_I y QR_es_A. Si las columnas son
    dependientes, ValueError."""
    gs = gram_schmidt(vectores)
    if gs["dependientes"]:
        j = gs["dependientes"][0]
        raise ValueError(f"las columnas son dependientes; QR con R invertible no existe (v{j} da u{j} = 0)")
    A = gs["A"]
    Q = sp.Matrix.hstack(*gs["q"])
    R = Q.T * A
    for i in range(R.rows):
        for j in range(i):
            if _es_cero(R[i, j]):
                R[i, j] = 0                      # abajo de la diagonal: cero exacto
    return {"Q": Q, "R": R, "A": A, "QtQ_es_I": columnas_ortonormales(Q),
            "QR_es_A": all(_es_cero(e) for e in Q * R - A)}


def compara_numpy_qr(vectores):
    """Compara la QR exacta con la de numpy (np.linalg.qr, reflexiones de Householder). Diccionario con Q_numpy,
    R_numpy, Q_exacta y R_exacta (en punto flotante), signos (el signo de cada entrada de la diagonal de R_numpy),
    coinciden (True si, al cambiar esos signos, numpy da la misma Q y la misma R, con tolerancia 1e-12 relativa)
    y texto (la explicación)."""
    ex = factorizacion_qr(vectores)
    Af = _flotante(ex["A"])
    Qn, Rn = np.linalg.qr(Af)
    Qe, Re = _flotante(ex["Q"]), _flotante(ex["R"])
    signos = np.where(np.diag(Rn) < 0, -1.0, 1.0)
    escala = max(1.0, float(np.abs(Af).max()))
    coinciden = bool(np.allclose(Qn * signos, Qe, rtol=0, atol=1e-12)
                     and np.allclose(signos[:, None] * Rn, Re, rtol=0, atol=1e-12 * escala))
    cambiadas = [j + 1 for j in range(len(signos)) if signos[j] < 0]
    if not cambiadas:
        texto = "numpy dio la misma factorización, con la diagonal de R positiva."
    else:
        cols = _lista(cambiadas)
        texto = (f"numpy usa reflexiones de Householder y aquí dio {'la columna' if len(cambiadas) == 1 else 'las columnas'} "
                 f"{cols} de Q con el signo cambiado, y {'el renglón' if len(cambiadas) == 1 else 'los renglones'} {cols} "
                 "de R también (por eso su diagonal tiene entradas negativas). Al cambiar esos signos, el producto QR no "
                 "cambia y queda la factorización exacta, con la diagonal de R positiva.")
    if not coinciden:
        texto += (" Aun con los signos ajustados quedan diferencias mayores que 1e-12: los vectores están cerca de ser "
                  "dependientes y el redondeo pesa.")
    return {"Q_numpy": Qn, "R_numpy": Rn, "Q_exacta": Qe, "R_exacta": Re, "signos": signos, "coinciden": coinciden,
            "texto": texto}


def _informe_gs(gs, normalizar, n_cifras):
    A = gs["A"]
    _nombra_vectores(A)
    pasos = gs["pasos"]
    k = A.cols
    for j in range(1, k + 1):
        cs = gs["coeficientes"][j]
        if cs:
            print("  " + ",  ".join(f"(v{j}·u{i})/(u{i}·u{i}) = {_valor(c, n_cifras)}" for i, c in cs))
        print(pasos[j - 1])
    if not gs["u"]:
        print("Todos los vectores son el cero: no queda ningún u." if k > 1 else "v1 es el vector cero: no queda ningún u.")
        return
    quedan = len(gs["u"])
    if gs["dependientes"]:
        dep = gs["dependientes"]
        print(f"Se {'descartó' if len(dep) == 1 else 'descartaron'} {_lista([f'v{j}' for j in dep])}: " +
              ("queda 1 vector, que genera" if quedan == 1 else f"quedan {quedan} vectores ortogonales, que generan") +
              f" lo mismo que los {k} dados.")
    nombres = [f"u{i}" for i in gs["indices"]]
    muestra(f"matriz de Gram de {_coma(nombres)}", gs["gram"], n_cifras)
    diagonal = all(gs["gram"][i, j] == 0 for i in range(len(nombres)) for j in range(len(nombres)) if i != j)
    if quedan == 1:
        print(f"Queda un solo vector, {nombres[0]}: no hay pares que revisar.")
    else:
        print(f"La matriz de Gram es diagonal: {_productos_nulos(nombres)}." if diagonal else
              "¡La matriz de Gram NO es diagonal!")
    if normalizar:
        for linea in pasos[k:]:
            print(linea)
        mide = all(_es_cero(_punto(q, q) - 1) for q in gs["q"])
        print((f"q{gs['indices'][0]} mide 1: " if quedan == 1 else "Cada q mide 1: ") + ("sí." if mide else "¡no!"))


def _informe_qr(gs, n_cifras):
    if gs["dependientes"]:
        j = gs["dependientes"][0]
        print(f"QR: las columnas de A son dependientes (v{j} dio u{j} = 0), así que no hay QR con R invertible.")
        return
    qr = factorizacion_qr(gs["A"])
    k = gs["A"].cols
    muestra("Q", qr["Q"], n_cifras)
    muestra("R = QᵀA", qr["R"], n_cifras, tex=r"R=Q^TA")
    diag = ", ".join(f"r{j}{j} = ‖u{j}‖" for j in range(1, k + 1)) if k < 10 else "cada r de la diagonal es la norma de su u"
    print(f"QᵀQ = I: {'sí' if qr['QtQ_es_I'] else 'no'};  QR = A: {'sí' if qr['QR_es_A'] else 'no'};  "
          + ("R es triangular superior con la diagonal positiva" if k > 1 else "R tiene una sola entrada, positiva")
          + f" ({diag}).")
    cmp = compara_numpy_qr(gs["A"])
    print("Comparación con np.linalg.qr (decimales): " + cmp["texto"])
    muestra_decimales("R de numpy", cmp["R_numpy"], n_cifras, tex=r"R_{\text{numpy}}")


def _gs_clasico_flotante(V):
    """Gram-Schmidt clásico en punto flotante (solo para la nota numérica): devuelve Q con columnas 'ortonormales'."""
    Q = []
    for v in V.T:
        w = v - sum((q @ v) * q for q in Q)
        Q.append(w / np.linalg.norm(w))
    return np.column_stack(Q)

# %% [markdown]
# **Recorrido paso a paso.** Con $\mathbf v_1=(1,1,0)$, $\mathbf v_2=(1,0,1)$ y $\mathbf v_3=(0,1,1)$: $\mathbf u_1=(1,1,0)$; $\mathbf u_2=\mathbf v_2-\tfrac12\mathbf u_1=(\tfrac12,-\tfrac12,1)$; $\mathbf u_3=\mathbf v_3-\tfrac12\mathbf u_1-\tfrac13\mathbf u_2=(-\tfrac23,\tfrac23,\tfrac23)$, y sus normalizados $\mathbf q_j$. Después, la QR de las dos primeras columnas, con $R=\begin{bmatrix}\sqrt2&\sqrt2/2\\0&\sqrt6/2\end{bmatrix}$, y la comparación con numpy. Al final, con $(2,1,1)$ en lugar de $(0,1,1)$, el tercer vector es $\mathbf v_1+\mathbf v_2$ y da $\mathbf u_3=\mathbf 0$.

# %%
_informe_gs(gram_schmidt("1, 1, 0; 1, 0, 1; 0, 1, 1", normalizar=True), True, 4)
print("\n" + "=" * 90 + "\n\nQR de las dos primeras columnas, A = [v1 v2]:")
_informe_qr(gram_schmidt("1, 1, 0; 1, 0, 1"), 4)
print("\n" + "=" * 90 + "\n\nCon (2, 1, 1) en lugar de (0, 1, 1):")
_informe_gs(gram_schmidt("1, 1, 0; 1, 0, 1; 2, 1, 1"), False, 4)

# %% [markdown]
# **La nota numérica, en números.** Los vectores $(1,\varepsilon,0,0)$, $(1,0,\varepsilon,0)$ y $(1,0,0,\varepsilon)$ con $\varepsilon=10^{-8}$ son independientes pero casi iguales. En punto flotante, $1+\varepsilon^2$ se redondea a $1$, y Gram-Schmidt clásico entrega $\mathbf q_2$ y $\mathbf q_3$ que forman un ángulo de unos $60^\circ$ en lugar de $90^\circ$. El cálculo exacto y `np.linalg.qr` no tienen ese problema.

# %%
_eps = sp.Rational(1, 10 ** 8)
_V = lista_vectores([(1, _eps, 0, 0), (1, 0, _eps, 0), (1, 0, 0, _eps)])
_Qc = _gs_clasico_flotante(_flotante(_V))
_Qn = np.linalg.qr(_flotante(_V))[0]
_qe = gram_schmidt(_V)["q"]
print(f"Gram-Schmidt clásico en decimales: q2·q3 ≈ {cifras(_Qc[:, 1] @ _Qc[:, 2], 4)}  (debería ser 0)")
print(f"np.linalg.qr (Householder):        |q2·q3| ≈ {cifras(abs(_Qn[:, 1] @ _Qn[:, 2]), 2)}")
print(f"Gram-Schmidt exacto (fracciones):  q2·q3 = {sp.simplify(_punto(_qe[1], _qe[2]))}")

# %% [markdown]
# **Calculadora.** Escribe los vectores **uno por renglón**, en el orden en que quieres ortogonalizarlos. La calculadora muestra cada paso, los $\mathbf u_j$, la matriz de Gram y, si lo pides, los $\mathbf q_j$ y la factorización $A=QR$ comparada con la de numpy. Prueba cambiar el orden de los vectores, agregar uno dependiente y usar vectores de $\mathbb R^4$.

# %%
def calculadora_gram_schmidt(vectores_txt, normalizar, mostrar_qr, n_cifras):
    with bloque():                 # el resultado aparece completo, de una sola vez
        try:
            gs = gram_schmidt(vectores_txt, normalizar, n_cifras)
        except ValueError as err:
            print("Revisa la entrada:", err); return
        _aviso_coma(vectores_txt)
        _informe_gs(gs, normalizar, n_cifras)
        if mostrar_qr:
            print()
            _informe_qr(gs, n_cifras)


_CONTROLES_3 = dict(
    vectores_txt=widgets.Textarea(value="1, 1, 0\n1, 0, 1\n0, 1, 1", description="vectores:", continuous_update=False,
                                  layout=widgets.Layout(width="340px", height="130px")),
    normalizar=widgets.Checkbox(value=True, description="normalizar (q)"),
    mostrar_qr=widgets.Checkbox(value=False, description="mostrar QR"),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_gram_schmidt, **_CONTROLES_3);

# %%
# Casos de prueba de la sección 3 (resultado conocido)
_GS3 = "1, 1, 0; 1, 0, 1; 0, 1, 1"             # los tres vectores del recorrido
PRUEBAS_3 = [
    ("gram_schmidt((1, 1, 0), (1, 0, 1)) → u1 = (1, 1, 0), u2 = (1/2, -1/2, 1)",
     lambda: (lambda g: g["u"][0] == vector("1,1,0") and g["u"][1] == vector("1/2, -1/2, 1"))(
         gram_schmidt("1,1,0; 1,0,1"))),
    ("con (0, 1, 1) como tercer vector: u3 = (-2/3, 2/3, 2/3) y la matriz de Gram es diag(2, 3/2, 4/3)",
     lambda: (lambda g: g["u"][2] == vector("-2/3, 2/3, 2/3")
                        and g["gram"] == sp.diag(2, sp.Rational(3, 2), sp.Rational(4, 3)))(gram_schmidt(_GS3))),
    ("con (2, 1, 1) como tercer vector: dependientes = [3] y quedan dos u",
     lambda: (lambda g: g["dependientes"] == [3] and len(g["u"]) == 2 and g["indices"] == [1, 2])(
         gram_schmidt("1,1,0; 1,0,1; 2,1,1"))),
    ("normalizar: q1 = (1/√2, 1/√2, 0) y cada q mide 1",
     lambda: (lambda g: g["q"][0] == sp.Matrix([1 / sp.sqrt(2), 1 / sp.sqrt(2), 0])
                        and all(sp.simplify(_punto(q, q)) == 1 for q in g["q"]))(gram_schmidt(_GS3, normalizar=True))),
    ("QR de [(1, 1, 0) (1, 0, 1)]: R = [[√2, √2/2], [0, √6/2]], QᵀQ = I y QR = A",
     lambda: (lambda f: sp.simplify(f["R"] - sp.Matrix([[sp.sqrt(2), sp.sqrt(2) / 2], [0, sp.sqrt(6) / 2]]))
                        == sp.zeros(2, 2) and f["QtQ_es_I"] and f["QR_es_A"]
                        and sp.simplify(f["Q"].T * f["Q"]) == sp.eye(2))(factorizacion_qr("1,1,0; 1,0,1"))),
    ("el orden importa: con (1, 0, 1), (1, 1, 0) sale u1 = (1, 0, 1) y otro u2, pero el mismo plano (rango 2)",
     lambda: (lambda g, h: g["u"][0] == vector("1,0,1") and g["u"][1] != h["u"][1]
                           and sp.Matrix.hstack(*g["u"], *h["u"]).rank() == 2)(
         gram_schmidt("1,0,1; 1,1,0"), gram_schmidt("1,1,0; 1,0,1"))),
    ("una base ya ortogonal, (1, 1, 1), (1, -1, 0), (1, 1, -2), no cambia",
     lambda: gram_schmidt("1,1,1; 1,-1,0; 1,1,-2")["u"] == [vector("1,1,1"), vector("1,-1,0"), vector("1,1,-2")]),
    ("en R⁴: gram_schmidt((1, 1, 0, 0), (1, 0, 1, 1)) → u2 = (1/2, -1/2, 1, 1)",
     lambda: gram_schmidt("1,1,0,0; 1,0,1,1")["u"][1] == vector("1/2, -1/2, 1, 1")),
    ("np.linalg.qr del recorrido: |Q de numpy| = |Q exacta| columna por columna (tolerancia 1e-12) y coinciden "
     "salvo signos",
     lambda: (lambda c: np.allclose(np.abs(c["Q_numpy"]), np.abs(c["Q_exacta"]), rtol=0, atol=1e-12)
                        and c["coinciden"])(compara_numpy_qr(_GS3))),
    ("gram_schmidt((0, 0, 0), (1, 0, 0)) → dependientes = [1] y u = [(1, 0, 0)]",
     lambda: (lambda g: g["dependientes"] == [1] and g["u"] == [vector("1,0,0")])(gram_schmidt("0,0,0; 1,0,0"))),
    ("factorizacion_qr((1, 2), (2, 4)) → ValueError «las columnas son dependientes»",
     lambda: _mensaje(lambda: factorizacion_qr("1,2; 2,4")).startswith("las columnas son dependientes")),
    ("prueba de humo: la calculadora corre con sus valores iniciales, también con QR y con un vector dependiente",
     lambda: _humo(calculadora_gram_schmidt, _CONTROLES_3)
             and _humo(calculadora_gram_schmidt, _CONTROLES_3, mostrar_qr=True)
             and _humo(calculadora_gram_schmidt, _CONTROLES_3, vectores_txt="1, 1, 0\n1, 0, 1\n2, 1, 1", mostrar_qr=True)),
    ("los textos de los pasos llevan el signo de la resta: «u2 = v2 - (1/2)·u1 = (1/2, -1/2, 1)»; con coeficiente "
     "negativo, (1, 0), (-1, 1) → «u2 = v2 + 1·u1 = (0, 1)»",
     lambda: gram_schmidt("1,1,0; 1,0,1")["pasos"][1] == "u2 = v2 - (1/2)·u1 = (1/2, -1/2, 1)"
             and gram_schmidt("1,0; -1,1")["pasos"][1] == "u2 = v2 + 1·u1 = (0, 1)"),
    ("compara_numpy_qr con columnas casi dependientes, (1, 1, 1), (1, 1, 1.00000001): coinciden = False",
     lambda: compara_numpy_qr("1,1,1; 1,1,1.00000001")["coinciden"] is False
             and compara_numpy_qr(_GS3)["coinciden"] is True),
    ("columnas_ortonormales: la Q del recorrido sí; [[1, 1], [0, 1], [0, 0]] y [[1, 0], [0, 2], [0, 0]] no",
     lambda: columnas_ortonormales(factorizacion_qr(_GS3)["Q"])
             and not columnas_ortonormales(sp.Matrix([[1, 1], [0, 1], [0, 0]]))
             and not columnas_ortonormales(sp.Matrix([[1, 0], [0, 2], [0, 0]]))),
]
corre_pruebas(PRUEBAS_3, 3)

# %% [markdown]
# **Contrasta (sección 3).** Aplica a mano Gram-Schmidt a $\mathbf v_1=(1,0,1)$ y $\mathbf v_2=(0,1,1)$, sin normalizar: $\mathbf u_1=\mathbf v_1$ y $\mathbf u_2=\mathbf v_2-\dfrac{\mathbf v_2\cdot\mathbf u_1}{\mathbf u_1\cdot\mathbf u_1}\,\mathbf u_1$. Escribe `mi_u2 = (a, b, c)`; puedes usar fracciones como texto (`"1/3"`) o decimales.

# %%
mi_u2 = None     # tres números, por ejemplo (1, -2, 0.5) o ("1/3", "2/3", 1)

if mi_u2 is None:     # la referencia se muestra después de tu cálculo
    print("Falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda.")
else:
    try:
        mio = vector(mi_u2, "tu u2")
    except ValueError as err:
        print(f"Revisa tu u2: {str(err).rstrip('.')}.")
    else:
        if mio.rows != 3:          # respuesta imposible: la referencia todavía no se muestra
            print(f"v1 y v2 son de R³, así que u2 tiene tres componentes (escribiste {mio.rows}): revisa tu respuesta "
                  "y vuelve a ejecutar la celda.")
        elif np.linalg.norm(_flotante(mio).ravel()) < 0.01 * np.linalg.norm(_flotante(vector("0, 1, 1")).ravel()):
            # imposible: v1 y v2 no son paralelos, así que u2 no es (casi) cero
            print("Tu u2 es (casi) el vector cero, y u2 = 0 solo pasa si v2 es múltiplo de v1, que aquí no lo es: "
                  "revisa tu cálculo y vuelve a ejecutar la celda.")
        else:
            ref_gs = gram_schmidt("1, 0, 1; 0, 1, 1")
            ref, v1 = ref_gs["u"][1], ref_gs["A"][:, 0]
            print(f"La calculadora da: u2 = {_tupla(ref, 4)}.")
            mf, rf, v1f = _flotante(mio).ravel(), _flotante(ref).ravel(), _flotante(v1).ravel()
            producto = float(mf @ v1f)
            factor = (mf @ rf) / (rf @ rf)
            multiplo = np.linalg.norm(mf - factor * rf) <= 0.01 * np.linalg.norm(mf)          # tolerancia relativa
            if np.all(np.abs(mf - rf) <= 0.01):
                print("Coinciden (con tolerancia de 0.01 por componente).")
            elif multiplo and abs(factor - 1) <= 0.02:
                print("Coinciden salvo redondeo (tu vector difiere de u2 en menos de 2 %; escribe las fracciones exactas "
                      "si puedes).")
            elif multiplo:
                print(f"Tu vector es un múltiplo de la respuesta del método ({cifras(factor, 4)} veces u2): también "
                      "sirve, porque es perpendicular a v1 y genera, con v1, el mismo plano. El método, sin escalar, "
                      "da el u2 de arriba.")
            elif abs(producto) > 0.01 * np.linalg.norm(mf) * np.linalg.norm(v1f):
                print(f"NO coinciden: tu u2 no es perpendicular a v1, porque u2·v1 ≈ {cifras(producto, 3)} (debería ser 0). "
                      "Revisa el coeficiente (v2·u1)/(u1·u1) y la resta.")
            else:
                print("NO coinciden: tu vector es perpendicular a v1, pero no está en el plano de v1 y v2. u2 debe ser "
                      "v2 menos un múltiplo de u1.")
            print("\nEl procedimiento completo:")
            _informe_gs(ref_gs, False, 4)
