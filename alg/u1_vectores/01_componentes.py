# ID: ALG-U1-NB01
# Notebook: alg/u1_vectores.ipynb · sección 1.1 vectores en 2D y 3D: componentes y magnitud
# Repositorio: alg/u1_vectores/01_componentes.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Vectores en 2D y 3D: componentes y magnitud
#
# El vector que va del punto $P$ al punto $Q$ se obtiene restando **llegada menos salida**:
#
# $$\overrightarrow{PQ}=Q-P=(q_x-p_x,\;q_y-p_y,\;q_z-p_z).$$
#
# Su magnitud (norma) es su longitud, por Pitágoras (en el espacio, Pitágoras se aplica dos veces: primero en el piso y luego en la vertical):
#
# $$\|\mathbf v\|=\sqrt{v_x^2+v_y^2}\quad\text{(plano)},\qquad \|\mathbf v\|=\sqrt{v_x^2+v_y^2+v_z^2}\quad\text{(espacio)}.$$
#
# En el plano, el ángulo $\theta$ que forma $\mathbf v$ con el eje $x$ cumple $\tan\theta=v_y/v_x$, pero despejar $\theta=\arctan(v_y/v_x)$ **pierde el cuadrante**: $(3,4)$ y $(-3,-4)$ dan el mismo cociente y apuntan en sentidos opuestos. La función de dos argumentos `arctan2(vy, vx)` mira el signo de cada componente y devuelve un ángulo en $(-180^\circ,\,180^\circ]$, o en $(-\pi,\,\pi]$ si trabajas en radianes. El vector cero no tiene dirección, así que no tiene ángulo.
#
# La calculadora recibe $P$ y $Q$ (2 o 3 componentes), escribe $\overrightarrow{PQ}$ y su magnitud y, en el plano, el ángulo con la unidad a la vista. Dibuja la flecha con sus componentes punteadas.

# %%
def componentes(P, Q):
    """Vector de P a Q (llegada menos salida). P y Q deben tener 2 o 3 componentes, las mismas los dos."""
    P, Q = vector(P, nombre="P"), vector(Q, nombre="Q")
    if len(P) != len(Q) or len(P) not in (2, 3):
        raise ValueError(f"P y Q deben tener la misma cantidad de componentes, 2 (plano) o 3 (espacio); "
                         f"tienen {len(P)} y {len(Q)}.")
    return Q - P


def magnitud(v):
    """Longitud del vector: raíz de la suma de los cuadrados de sus componentes."""
    return _norma(vector(v))


def angulo_2d(v, unidad="grados"):
    """Ángulo de un vector del plano con el eje x, con arctan2: en (-180°, 180°] o (-π, π]."""
    v = vector(v)
    if len(v) != 2:
        raise ValueError("el ángulo con el eje x se calcula aquí solo para vectores del plano (2 componentes).")
    if not np.any(v):
        raise ValueError("el vector cero no tiene dirección.")
    t = math.atan2(v[1], v[0])
    if t == -math.pi:                         # (-1, -0.0) daría -180°: se reporta como 180°
        t = math.pi
    return en_unidad(t, unidad)


def calculadora_componentes(P_txt, Q_txt, unidad, n_cifras):
    try:
        P, Q = vector(P_txt, nombre="P"), vector(Q_txt, nombre="Q")
        v = componentes(P, Q)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    cuadrados = " + ".join(f"{c}²" for c in ("vx", "vy", "vz")[:len(v)])
    print(f"P = {fmt(P, n_cifras)}    Q = {fmt(Q, n_cifras)}")
    print(f"PQ = Q − P = {fmt(v, n_cifras)}      (llegada menos salida)")
    print(f"‖PQ‖ = √({cuadrados}) = √{cifras(math.fsum(v**2), n_cifras)} = {cifras(magnitud(v), n_cifras)}")
    if not np.any(v):
        print("P = Q: PQ es el vector cero, que no tiene dirección (ni ángulo ni vector unitario).")
    elif len(v) == 2:
        try:
            ang = angulo_2d(v, unidad)
            print(f"ángulo con el eje x (arctan2) = {cifras(ang, n_cifras)}{_sufijo(unidad)}")
            if v[0] < 0:
                malo = en_unidad(math.atan(v[1] / v[0]), unidad)
                print(f"Ojo: arctan(vy/vx) daría {cifras(malo, n_cifras)}{_sufijo(unidad)}, el ángulo del vector "
                      "opuesto: arctan pierde el cuadrante.")
        except ValueError as err:
            print("ángulo:", err)
    else:
        print("En 3D un solo ángulo no basta para fijar la dirección; la da el vector unitario (sección 2).")
    with np.errstate(all="ignore"):          # componentes cercanas a 1e100: matplotlib eleva al cuadrado
        _dibuja_pq(P, Q, v)


def _dibuja_pq(P, Q, v):
    fig = plt.figure(figsize=(5.5, 4.5))
    if len(v) == 2:
        ax = fig.add_subplot()
        ax.plot([P[0], Q[0]], [P[1], P[1]], ":", color="tab:red", lw=1.5, label=f"componente x = {cifras(v[0], 3)}")
        ax.plot([Q[0], Q[0]], [P[1], Q[1]], ":", color="tab:green", lw=1.5, label=f"componente y = {cifras(v[1], 3)}")
        ax.quiver(*P, *v, angles="xy", scale_units="xy", scale=1, color="navy", width=0.012)
        ax.legend(loc="best", fontsize=8); ax.grid(True)
    else:
        ax = fig.add_subplot(projection="3d")
        esquinas = np.array([P, (Q[0], P[1], P[2]), (Q[0], Q[1], P[2]), Q])
        ax.plot(*esquinas.T, "k:", lw=1.2, label="componentes x, y, z")
        ax.quiver(*P, *v, color="navy", arrow_length_ratio=0.12, lw=2)
        ax.set_zlabel("z"); ax.legend(loc="upper left", fontsize=8)
        _vista(ax, [v])                        # que la flecha no apunte hacia el ojo
    for X, e in (((P, "P = Q"),) if not np.any(v) else ((P, "P"), (Q, "Q"))):
        ax.scatter(*X, color="black", s=15)
        ax.text(*X, f"  {e}", fontsize=11)
    _ejes_iguales(ax, [P, Q])
    ax.set_xlabel("x"); ax.set_ylabel("y")
    ax.set_title(f"PQ = {fmt(v, 3)},  ‖PQ‖ = {cifras(magnitud(v), 3)}")
    plt.show()


widgets.interact(calculadora_componentes,
    P_txt=widgets.Text(value="1, 1", description="P =", continuous_update=False),
    Q_txt=widgets.Text(value="4, 5", description="Q =", continuous_update=False),
    unidad=widgets.Dropdown(options=["grados", "radianes"], value="grados", description="unidad"),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras sig.", continuous_update=False));

# %%
# Casos de prueba de la sección 1 (resultado conocido)
PRUEBAS_1 = [
    ("PQ de (1, 1) a (4, 5) = (3, 4), magnitud 5, ángulo 53.13°",
     lambda: np.allclose(componentes((1, 1), (4, 5)), (3, 4)) and cerca(magnitud(componentes((1, 1), (4, 5))), 5)
             and round(angulo_2d(componentes((1, 1), (4, 5)), "grados"), 2) == 53.13),
    ("PQ de (1, -2, 4) a (3, 4, 1) = (2, 6, -3), magnitud 7",
     lambda: np.allclose(componentes((1, -2, 4), (3, 4, 1)), (2, 6, -3))
             and cerca(magnitud(componentes((1, -2, 4), (3, 4, 1))), 7)),
    ("ángulo de (-3, -4) = -126.87° (no 53.13°)", lambda: round(angulo_2d((-3, -4), "grados"), 2) == -126.87),
    ("el vector (0, 0) no tiene ángulo", lambda: _rechaza(lambda: angulo_2d((0, 0)))),
    ("entradas mal escritas se rechazan: '1, a' y vacía",
     lambda: _rechaza(lambda: vector("1, a")) and _rechaza(lambda: vector(""))),
]
for nombre, prueba in PRUEBAS_1:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 1).** Calcula a mano $\|\overrightarrow{PQ}\|$ para $P(-1,2,3)$ y $Q(1,-1,9)$: primero las componentes (llegada menos salida), después la raíz de la suma de los cuadrados. Escribe tu resultado.

# %%
mi_magnitud = None     # escribe tu resultado, por ejemplo: 5.385 o "sqrt(29)"

if mi_magnitud is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    v_ref = componentes((-1, 2, 3), (1, -1, 9))
    print("La calculadora da: PQ =", fmt(v_ref, 4), "  ‖PQ‖ =", cifras(magnitud(v_ref), 4))
    try:
        ok = cerca(numero(mi_magnitud, "tu resultado"), magnitud(v_ref), 1e-3)
    except ValueError as err:
        print("Revisa tu valor:", err); ok = None
    if ok is not None:
        print("coinciden" if ok else
              "NO coinciden: revisa las restas (llegada menos salida) y eleva al cuadrado también las componentes negativas")
