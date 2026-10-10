# ID: ALG-U1-NB04
# Notebook: alg/u1_vectores.ipynb · sección 1.4 producto cruz
# Repositorio: alg/u1_vectores/04_cruz.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 4. Producto cruz
#
# El producto cruz combina dos vectores **de 3 componentes** en un tercer vector (definición 1.13 del libro):
#
# $$\mathbf a\times\mathbf b=\big(a_yb_z-a_zb_y,\;\;a_zb_x-a_xb_z,\;\;a_xb_y-a_yb_x\big).$$
#
# * Es **perpendicular** a los dos factores: $(\mathbf a\times\mathbf b)\cdot\mathbf a=0$ y $(\mathbf a\times\mathbf b)\cdot\mathbf b=0$. Su sentido lo da la regla de la mano derecha.
# * Su magnitud es el **área** del paralelogramo que forman: $\|\mathbf a\times\mathbf b\|=\|\mathbf a\|\,\|\mathbf b\|\,\mathrm{sen}\,\theta$. La mitad es el área del triángulo.
# * **El orden importa:** $\mathbf b\times\mathbf a=-(\mathbf a\times\mathbf b)$ (anticonmutatividad). Un momento $\mathbf r\times\mathbf F$ escrito al revés apunta en sentido contrario.
# * $\mathbf a\times\mathbf b=\mathbf 0$ si los vectores son **paralelos** (o alguno es el vector cero): no encierran área.
# * Solo existe en 3D. Si das dos vectores del plano, la calculadora les agrega $z=0$ y te avisa; el resultado queda sobre el eje $z$.

# %%
def _cruz(a, b):
    """a × b de dos vectores de 3 componentes ya leídos, con la cancelación por redondeo limpiada."""
    terminos = np.array([[a[1] * b[2], a[2] * b[1]],
                         [a[2] * b[0], a[0] * b[2]],
                         [a[0] * b[1], a[1] * b[0]]])
    return _limpia(terminos[:, 0] - terminos[:, 1], np.abs(terminos).sum(axis=1))


def cruz(a, b):
    """(a × b, aviso). aviso es '' en 3D o un mensaje si a vectores del plano se les agregó z = 0."""
    a, b = _par(a, b)
    aviso = ""
    if len(a) == 2:
        a, b = np.append(a, 0.0), np.append(b, 0.0)
        aviso = "Aviso: el producto cruz solo existe en 3D; a los vectores del plano se les agregó z = 0."
    elif len(a) != 3:
        raise ValueError(f"el producto cruz solo está definido para vectores de 3 componentes "
                         f"(o del plano, agregando z = 0); estos tienen {len(a)}.")
    return _cruz(a, b), aviso


def paralelos(a, b, tol=1e-12):
    """True si a × b = 0 dentro del redondeo: mismo sentido, sentido opuesto, o alguno es el vector cero."""
    a, b = _par(a, b)
    if not np.any(a) or not np.any(b):
        return True
    ea, eb = _escalado(a), _escalado(b)       # escalar por potencias de 2 no cambia la dirección
    return _norma(cruz(ea, eb)[0]) <= tol * _norma(ea) * _norma(eb)


def calculadora_cruz(a_txt, b_txt, n_cifras):
    try:
        a, b = _par(a_txt, b_txt)
        c, aviso = cruz(a, b)
        d, _ = cruz(b, a)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    if aviso:
        print(aviso)
        a, b = np.append(a, 0.0), np.append(b, 0.0)
    area = _norma(c)
    print(f"a × b = {fmt(a, n_cifras)} × {fmt(b, n_cifras)} = {fmt(c, n_cifras)}")
    print(f"b × a = {fmt(b, n_cifras)} × {fmt(a, n_cifras)} = {fmt(d, n_cifras)}   (el orden importa: b × a = −(a × b))")
    print(f"área del paralelogramo = ‖a × b‖ = {cifras(area, n_cifras)}")
    print(f"área del triángulo = ‖a × b‖ / 2 = {cifras(area / 2, n_cifras)}")
    print(f"comprobación: (a × b) · a = {cifras(punto(c, a), n_cifras)}   y   (a × b) · b = {cifras(punto(c, b), n_cifras)}",
          end="")
    if paralelos(a, b):
        print("\na × b = 0: los vectores son paralelos (o alguno es el vector cero) y no encierran área.")
    else:
        print("   (a × b es perpendicular a los dos)")
    with np.errstate(all="ignore"):          # componentes cercanas a 1e100: matplotlib eleva al cuadrado
        _dibuja_cruz(a, b, c, area)


def _dibuja_cruz(a, b, c, area):
    o = np.zeros(3)
    fig = plt.figure(figsize=(6, 5.5))
    ax = fig.add_subplot(projection="3d")
    for v, color, e in ((a, "tab:red", "a"), (b, "tab:green", "b"), (c, "navy", "a × b")):
        ax.quiver(*o, *v, color=color, arrow_length_ratio=0.1, lw=2, label=e)
    ax.plot(*np.array([a, a + b, b]).T, "k:", lw=1, label="paralelogramo")
    _ejes_iguales(ax, [o, a, b, a + b, c])
    _vista(ax, [a, b], normal=c)               # ninguna flecha apunta al ojo; el paralelogramo, inclinado
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.set_title(f"a × b = {fmt(c, 3)}, área = {cifras(area, 3)}")
    ax.legend(loc="upper left", fontsize=8)
    plt.show()


widgets.interact(calculadora_cruz,
    a_txt=widgets.Text(value="1, 2, 3", description="a =", continuous_update=False),
    b_txt=widgets.Text(value="4, 5, 6", description="b =", continuous_update=False),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras sig.", continuous_update=False));

# %%
# Casos de prueba de la sección 4 (resultado conocido)
PRUEBAS_4 = [
    ("(1, 2, 3) × (4, 5, 6) = (-3, 6, -3), en orden invertido (3, -6, 3); no son paralelos",
     lambda: np.allclose(cruz((1, 2, 3), (4, 5, 6))[0], (-3, 6, -3)) and np.allclose(cruz((4, 5, 6), (1, 2, 3))[0], (3, -6, 3))
             and not paralelos((1, 2, 3), (4, 5, 6))),
    ("i × j = k y j × i = -k",
     lambda: np.allclose(cruz((1, 0, 0), (0, 1, 0))[0], (0, 0, 1)) and np.allclose(cruz((0, 1, 0), (1, 0, 0))[0], (0, 0, -1))),
    ("(1, -3, 2) × (-2, 6, -4) = 0: son paralelos",
     lambda: np.allclose(cruz((1, -3, 2), (-2, 6, -4))[0], 0) and paralelos((1, -3, 2), (-2, 6, -4))),
    ("(3, 0) × (0, 4) = (0, 0, 12), con aviso de z = 0",
     lambda: np.allclose(cruz((3, 0), (0, 4))[0], (0, 0, 12)) and cruz((3, 0), (0, 4))[1] != ""),
    ("vectores de 4 componentes se rechazan", lambda: _rechaza(lambda: cruz((1, 2, 3, 4), (5, 6, 7, 8)))),
]
for nombre, prueba in PRUEBAS_4:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 4).** Calcula a mano $(2,-1,0)\times(1,3,2)$ con la fórmula de la definición 1.13 y escribe el resultado. Después invierte el orden y calcula $(1,3,2)\times(2,-1,0)$.

# %%
mi_axb = None     # tu a × b como texto, por ejemplo: "1, -2, 3"
mi_bxa = None     # después, b × a (el orden invertido)

a_c, b_c = (2, -1, 0), (1, 3, 2)
if mi_axb is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    print("La calculadora da: a × b =", fmt(cruz(a_c, b_c)[0], 4))
    try:
        ok = np.allclose(vector(mi_axb, 3, "tu a × b"), cruz(a_c, b_c)[0], rtol=0, atol=1e-6)
    except ValueError as err:
        print("Revisa tu valor:", err); ok = None
    if ok is not None:
        print("a × b: coinciden" if ok else
              "a × b: NO coinciden; usa (ay bz − az by, az bx − ax bz, ax by − ay bx) y cuida los signos")
    if mi_bxa is None:
        print("falta tu b × a: escríbelo arriba y vuelve a ejecutar la celda")
    else:
        print("La calculadora da: b × a =", fmt(cruz(b_c, a_c)[0], 4))
        try:
            ok = np.allclose(vector(mi_bxa, 3, "tu b × a"), cruz(b_c, a_c)[0], rtol=0, atol=1e-6)
        except ValueError as err:
            print("Revisa tu valor:", err); ok = None
        if ok is not None:
            print("b × a: coinciden; es −(a × b)" if ok else
                  "b × a: NO coinciden; al invertir el orden, cada componente cambia de signo")
