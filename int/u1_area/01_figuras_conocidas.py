# ID: INT-U1-NB01
# Notebook: int/u1_area.ipynb · sección 1.1 áreas de figuras conocidas
# Repositorio: int/u1_area/01_figuras_conocidas.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Áreas de figuras conocidas: área bajo una poligonal
#
# Si la gráfica de $f\geq0$ está formada por segmentos de recta, la región bajo ella se parte en **trapecios**, uno por cada par de puntos consecutivos $(x_k,y_k)$ y $(x_{k+1},y_{k+1})$:
#
# $$A_k=\frac{y_k+y_{k+1}}{2}\,(x_{k+1}-x_k).$$
#
# Un rectángulo es un trapecio con lados paralelos iguales y un triángulo, uno con un lado paralelo de longitud cero. Por **aditividad**, el área total es la suma de las $A_k$. Las unidades del área son (unidad del eje $y$)·(unidad del eje $x$).

# %%
def area_poligonal(xs, ys):
    """Áreas de los trapecios bajo la poligonal que une (xs[k], ys[k]) y su suma.
    Exige al menos dos puntos, x estrictamente creciente y y >= 0."""
    xs, ys = np.asarray(xs, dtype=float), np.asarray(ys, dtype=float)
    if xs.size != ys.size:
        raise ValueError(f"hay {xs.size} valores de x y {ys.size} de y: deben ser los mismos.")
    if xs.size < 2:
        raise ValueError("se necesitan al menos dos puntos para formar un trapecio.")
    if np.any(np.diff(xs) <= 0):
        raise ValueError("los valores de x deben ir de menor a mayor, sin repetirse.")
    if np.any(ys < 0):
        raise ValueError("hay valores de y negativos: esta sección mide el área de regiones sobre el eje. "
                         "Para áreas con signo usa la calculadora de la sección 3.")
    areas = (ys[:-1] + ys[1:]) / 2 * np.diff(xs)
    return areas, float(areas.sum())


def calculadora_poligonal(x_txt, y_txt, unidad_x, unidad_y, n_cifras):
    try:
        xs, ys = lista_numeros(x_txt, "los valores de x"), lista_numeros(y_txt, "los valores de y")
        areas, total = area_poligonal(xs, ys)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    if unidad_x and unidad_x == unidad_y:
        u = f" {unidad_x}²"
    elif unidad_x or unidad_y:
        u = " " + "·".join(v for v in (unidad_y, unidad_x) if v)
    else:
        u = ""
    for k, a in enumerate(areas):
        print(f"  trapecio {k + 1}: x de {cifras(xs[k], n_cifras)} a {cifras(xs[k + 1], n_cifras)} -> área {cifras(a, n_cifras)}{u}")
    print(f"área total ≈ {cifras(total, n_cifras)}{u}   ({n_cifras} cifras significativas, redondeado)")
    fig, ax = plt.subplots(figsize=(5, 3))
    ax.fill_between(xs, 0, ys, color="0.85")
    ax.plot(xs, ys, "o-", color="navy")
    for xv in xs[1:-1]:
        ax.axvline(xv, color="gray", ls="--", lw=0.7)
    ax.axhline(0, color="black", lw=0.7)
    ax.set_xlabel(f"x ({unidad_x})" if unidad_x else "x"); ax.set_ylabel(f"y ({unidad_y})" if unidad_y else "y")
    ax.set_title("Área bajo la poligonal"); plt.show()


widgets.interact(calculadora_poligonal,
    x_txt=widgets.Text(value="0, 1, 3, 4", description="x ="),
    y_txt=widgets.Text(value="1, 3, 3, 0", description="y ="),
    unidad_x=widgets.Text(value="", description="unidad x"),
    unidad_y=widgets.Text(value="", description="unidad y"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 1 (resultado conocido)
PRUEBAS_1 = [
    ("rectángulo: f = 3 en [1, 5] da 12", lambda: cerca(area_poligonal([1, 5], [3, 3])[1], 12)),
    ("trapecio: f = 2x + 1 en [0, 3] da 12", lambda: cerca(area_poligonal([0, 3], [1, 7])[1], 12)),
    ("poligonal de la figura 1.1(c): 2 + 6 + 1.5 = 9.5",
     lambda: np.allclose(area_poligonal([0, 1, 3, 4], [1, 3, 3, 0])[0], [2, 6, 1.5])),
    ("|x - 2| en [0, 5]: dos triángulos, 6.5", lambda: cerca(area_poligonal([0, 2, 5], [2, 0, 3])[1], 6.5)),
    ("x decreciente o y negativa se rechazan", lambda: _rechaza(lambda: area_poligonal([0, 2, 1], [1, 1, 1]))
                                                     and _rechaza(lambda: area_poligonal([0, 1], [1, -1]))),
]
for nombre, prueba in PRUEBAS_1:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 1).** Calcula a mano el área bajo la poligonal que une $(0,0)$, $(2,4)$, $(5,4)$ y $(6,0)$. Escribe tu resultado y ejecuta la celda.

# %%
mi_area = None           # escribe un número, por ejemplo: 12.5

if mi_area is None:      # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    ref = area_poligonal([0, 2, 5, 6], [0, 4, 4, 0])[1]
    print("La calculadora da:", cifras(ref, 4))
    print("coinciden" if cerca(mi_area, ref, 1e-6) else
          "NO coinciden: parte la región en un triángulo, un rectángulo y otro triángulo")
