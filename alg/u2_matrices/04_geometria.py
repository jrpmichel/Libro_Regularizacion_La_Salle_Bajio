# ID: ALG-U2-NB04
# Notebook: alg/u2_matrices.ipynb · sección 2.4 interpretación geométrica del determinante
# Repositorio: alg/u2_matrices/04_geometria.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 4. Interpretación geométrica del determinante
#
# Una matriz $A$ de $2\times 2$ transforma cada punto $\mathbf x$ del plano en $A\mathbf x$. Las **columnas** de $A$ son las imágenes de los vectores de la base:
#
# $$A=\begin{bmatrix}a_{11}&a_{12}\\a_{21}&a_{22}\end{bmatrix},\qquad A\hat{\mathbf i}=\begin{bmatrix}a_{11}\\a_{21}\end{bmatrix},\qquad A\hat{\mathbf j}=\begin{bmatrix}a_{12}\\a_{22}\end{bmatrix}.$$
#
# - El cuadrado unitario (lados $\hat{\mathbf i}$ y $\hat{\mathbf j}$) se transforma en el paralelogramo de lados $A\hat{\mathbf i}$ y $A\hat{\mathbf j}$, de área $|\det A|$. Cualquier figura cambia su área por el mismo factor $|\det A|$.
# - Si $\det A>0$, la orientación se **conserva**: lo que se recorría en sentido antihorario sigue en sentido antihorario.
# - Si $\det A<0$, la orientación se **invierte**: la figura queda reflejada, como en un espejo (una letra F se ve al revés).
# - Si $\det A=0$, el plano se **aplasta** sobre una recta (o sobre un punto, si $A$ es la matriz cero): el área de la imagen es 0 y la transformación no se puede deshacer, por eso $A$ no tiene inversa.
#
# En el espacio, el volumen del paralelepípedo de aristas $\mathbf a$, $\mathbf b$, $\mathbf c$ es el valor absoluto del determinante que las tiene como columnas, que es el triple producto escalar:
#
# $$V=\left|\det\,[\,\mathbf a\ \ \mathbf b\ \ \mathbf c\,]\right|=\left|\mathbf a\cdot(\mathbf b\times\mathbf c)\right|.$$
#
# Mueve las entradas de $A$ con los controles. La figura original aparece punteada y su imagen rellena; un triángulo hueco marca el punto $(1,0)$, la punta de $\hat{\mathbf i}$, y el triángulo naranja marca su imagen, la punta de $A\hat{\mathbf i}$, para que veas si la figura dio la vuelta.

# %%
_FIGURAS = {
    "cuadrado unitario": [(0, 0), (1, 0), (1, 1), (0, 1)],
    "letra F": [(0, 0), (0.2, 0), (0.2, 0.45), (0.5, 0.45), (0.5, 0.6), (0.2, 0.6), (0.2, 0.85), (0.6, 0.85), (0.6, 1), (0, 1)],
}
_AZUL, _NARANJA, _VERDE = "#2a78d6", "#eb6834", "#1baf7a"


def efecto(A):
    """Efecto de una matriz 2×2 sobre el plano: det, factor de área |det|, orientación e imagen del plano."""
    A = como_matriz(A, "A")
    if A.shape != (2, 2):
        raise ValueError(f"esta sección usa matrices 2×2 (esta es {dims(A)}); para tres vectores usa volumen(a, b, c)")
    d = A.det()
    if d > 0:
        orientacion = "conserva"
    elif d < 0:
        orientacion = "invierte"
    else:
        orientacion = "aplasta"
    imagen_plano = "plano" if d != 0 else ("recta" if A.rank() == 1 else "punto")
    return {"det": d, "area": abs(d), "orientacion": orientacion, "imagen": imagen_plano}


def volumen(a, b, c):
    """|a·(b×c)|: volumen exacto del paralelepípedo de aristas a, b, c."""
    vs = [vector(v, nombre) for v, nombre in ((a, "a"), (b, "b"), (c, "c"))]
    for v, nombre in zip(vs, "abc"):
        if v.rows != 3:
            raise ValueError(f"{nombre} debe tener 3 componentes (tiene {v.rows})")
    a, b, c = vs
    return abs(a.dot(b.cross(c)))


def imagen(A, puntos):
    """Imagen A·p de cada punto p = (x, y), como arreglo de numpy (para graficar)."""
    M = np.array(como_matriz(A, "A").tolist(), dtype=float)
    return np.asarray(puntos, dtype=float) @ M.T


def area_poligono(puntos):
    """Área con signo (fórmula del cordón de zapato): positiva si los vértices van en sentido antihorario."""
    P = np.asarray(puntos, dtype=float)
    xs, ys = P[:, 0], P[:, 1]
    return 0.5 * float(np.dot(xs, np.roll(ys, -1)) - np.dot(np.roll(xs, -1), ys))


def _corto(v, parentesis=False):
    """Número exacto con a lo más dos decimales (los controles avanzan de 0.1 en 0.1); negativos entre paréntesis si se pide."""
    texto = f"{float(v):.2f}".rstrip("0").rstrip(".") if v != 0 else "0"
    return f"({texto})" if parentesis and v < 0 else texto


def calculadora_geometria(a11, a12, a21, a22, figura):
    A = como_matriz([[round(a11, 1), round(a12, 1)], [round(a21, 1), round(a22, 1)]], "A")
    e = efecto(A)
    P = np.array(_FIGURAS[figura], dtype=float)
    Q = imagen(A, P)
    i_img, j_img = imagen(A, [(1, 0), (0, 1)])
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    cerrar = lambda R: np.vstack([R, R[:1]])
    ax.plot(*cerrar(P).T, "--", color="0.45", lw=1.5, label=f"{figura} (original)")
    if e["imagen"] == "plano":
        ax.fill(*Q.T, color=_AZUL, alpha=0.3, lw=0)
        ax.plot(*cerrar(Q).T, color=_AZUL, lw=2, label="imagen bajo A")
    elif e["imagen"] == "recta":
        ax.plot(*cerrar(Q).T, color=_AZUL, lw=3, label="imagen: un segmento")
    else:
        ax.plot(0, 0, "o", color=_AZUL, ms=9, label="imagen: un punto")
    for v, color, texto in ((i_img, _NARANJA, "Aî"), (j_img, _VERDE, "Aĵ")):
        if np.hypot(*v) > 1e-12:
            ax.annotate("", xy=v, xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=color, lw=1.8))
            ax.annotate(texto, xy=v, xytext=(6, 6), textcoords="offset points", fontsize=10)
    ax.plot(1, 0, "^", ms=10, mfc="none", mec="0.45", ls="none", label="punto (1, 0) = punta de î")
    ax.plot(*i_img, "^", ms=11, color=_NARANJA, ls="none", label="imagen de (1, 0) = punta de Aî")
    lim = 1.15 * max(1.0, float(np.abs(np.vstack([P, Q, [i_img, j_img]])).max()))
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim); ax.set_aspect("equal")
    ax.axhline(0, color="black", lw=0.6); ax.axvline(0, color="black", lw=0.6)
    ax.grid(alpha=0.3); ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title(f"det A = {_corto(e['det'])}")
    ax.legend(loc="best", fontsize=8)
    plt.show()
    c = [[_corto(A[i, j]) for j in range(2)] for i in range(2)]
    p = [[_corto(A[i, j], parentesis=True) for j in range(2)] for i in range(2)]
    print(f"A = [{c[0][0]}, {c[0][1]}; {c[1][0]}, {c[1][1]}]   (columnas: Aî = ({c[0][0]}, {c[1][0]}), Aĵ = ({c[0][1]}, {c[1][1]}))")
    print(f"det A = {p[0][0]}·{p[1][1]} - {p[0][1]}·{p[1][0]} = {_corto(e['det'])}")
    if e["orientacion"] == "aplasta":
        if e["imagen"] == "recta":
            print("det A = 0: el plano se aplasta sobre una recta y la imagen de la figura es un segmento (área 0).")
        else:
            print("det A = 0 y A es la matriz cero: todo el plano cae en el origen; la imagen es un punto (área 0).")
        print("Como dos puntos distintos llegan al mismo lugar, A no se puede deshacer: no tiene inversa.")
        return
    print(f"Factor de área: |det A| = {_corto(e['area'])}. Área de la figura: {cifras(abs(area_poligono(P)), 4)} → "
          f"área de su imagen: {cifras(abs(area_poligono(Q)), 4)} (= {_corto(e['area'])} × {cifras(abs(area_poligono(P)), 4)}).")
    if e["orientacion"] == "conserva":
        print("Orientación: se conserva (det A > 0); el recorrido antihorario de la figura sigue siendo antihorario.")
    else:
        print("Orientación: se invierte (det A < 0); la figura queda reflejada, como vista en un espejo.")


_control = lambda valor, nombre: widgets.FloatSlider(value=valor, min=-3, max=3, step=0.1, description=nombre,
                                                     continuous_update=False, readout_format=".1f")
widgets.interact(calculadora_geometria,
    a11=_control(2, "a11"), a12=_control(1, "a12"), a21=_control(1, "a21"), a22=_control(1, "a22"),
    figura=widgets.Dropdown(options=list(_FIGURAS), value="cuadrado unitario", description="figura"));

# %%
# Casos de prueba de la sección 4 (resultado conocido)
_G = matriz("-1, 0.5; 0.3, 2")
PRUEBAS_4 = [
    ("[[2,1],[1,1]]: det 1, conserva", lambda: (lambda e: e["det"] == 1 and e["orientacion"] == "conserva")(efecto(matriz("2, 1; 1, 1")))),
    ("[[1,2],[2,4]]: det 0, aplasta sobre una recta",
     lambda: (lambda e: e["det"] == 0 and e["orientacion"] == "aplasta" and e["imagen"] == "recta")(efecto(matriz("1, 2; 2, 4")))),
    ("[[-1,0],[0,1]]: det -1, invierte", lambda: (lambda e: e["det"] == -1 and e["orientacion"] == "invierte")(efecto(matriz("-1, 0; 0, 1")))),
    ("[[1.2,0.4],[0,0.8]]: factor de área 0.96 exacto", lambda: efecto(matriz("1.2, 0.4; 0, 0.8"))["area"] == sp.Rational(24, 25)),
    ("matriz cero: la imagen es un punto", lambda: efecto(matriz("0, 0; 0, 0"))["imagen"] == "punto"),
    ("volumen((2,0,0),(0,3,0),(0,0,4)) = 24", lambda: volumen("2, 0, 0", "0, 3, 0", "0, 0, 4") == 24),
    ("volumen((1,0,0),(1,1,0),(1,1,1)) = 1", lambda: volumen((1, 0, 0), (1, 1, 0), (1, 1, 1)) == 1),
    ("área de la F = 0.305 y área con signo de su imagen = det A × 0.305",
     lambda: cerca(area_poligono(_FIGURAS["letra F"]), 0.305)
             and cerca(area_poligono(imagen(_G, _FIGURAS["letra F"])), float(efecto(_G)["det"]) * 0.305)),
    ("efecto rechaza una 3×3; volumen rechaza vectores de 2 componentes",
     lambda: _rechaza(lambda: efecto(matriz("1, 0, 0; 0, 1, 0; 0, 0, 1"))) and _rechaza(lambda: volumen("1, 0", "0, 1", "1, 1"))),
]
for nombre, prueba in PRUEBAS_4:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 4).** Calcula a mano el área de la imagen del cuadrado unitario bajo $A=\begin{bmatrix}3&1\\1&2\end{bmatrix}$ y escribe el número.

# %%
mi_area = None     # escribe un número (también sirve como texto, por ejemplo "2.5")

if mi_area is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    try:
        mio = numero(mi_area, "tu área")
    except ValueError as err:
        print("Revisa tu número:", err)
    else:
        ref = efecto(matriz("3, 1; 1, 2"))["area"]
        print("La calculadora da: área =", ref)
        print("coinciden" if mio == ref else
              "NO coinciden: el cuadrado unitario tiene área 1, así que su imagen tiene área |det A| = |a₁₁a₂₂ - a₁₂a₂₁|")
