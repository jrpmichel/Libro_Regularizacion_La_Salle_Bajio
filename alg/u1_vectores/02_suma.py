# ID: ALG-U1-NB02
# Notebook: alg/u1_vectores.ipynb · sección 1.2 suma, resta y vector unitario
# Repositorio: alg/u1_vectores/02_suma.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Suma, resta y vector unitario
#
# * **Suma** componente a componente: $\mathbf a+\mathbf b=(a_x+b_x,\;a_y+b_y,\;a_z+b_z)$. Solo se suman vectores con la misma cantidad de componentes. La **resta** es $\mathbf a-\mathbf b=\mathbf a+(-1)\,\mathbf b$.
# * **Múltiplo escalar:** $k\mathbf a=(ka_x,\;ka_y,\;ka_z)$ y $\|k\mathbf a\|=|k|\,\|\mathbf a\|$. Si $k<0$, el sentido se invierte.
# * **Vector unitario:** $\hat{\mathbf u}=\dfrac{\mathbf v}{\|\mathbf v\|}$ para $\mathbf v\neq\mathbf 0$. Tiene la misma dirección y sentido que $\mathbf v$ y magnitud 1. El vector cero no tiene vector unitario.
# * **Desigualdad del triángulo:** $\|\mathbf a+\mathbf b\|\le\|\mathbf a\|+\|\mathbf b\|$. La igualdad vale solo si los dos apuntan en el mismo sentido (o uno es cero). Por eso dos fuerzas perpendiculares de 300 N y 400 N dan una resultante de 500 N, no de 700 N.
#
# La calculadora forma $\alpha\mathbf a+\beta\mathbf b$ (con $\alpha=1$, $\beta=-1$ obtienes la resta), da su magnitud y su vector unitario, y compara $\|\alpha\mathbf a+\beta\mathbf b\|$ con $|\alpha|\,\|\mathbf a\|+|\beta|\,\|\mathbf b\|$. En el plano dibuja el paralelogramo.

# %%
def combinacion(a, b, alfa=1, beta=1):
    """α a + β b, componente a componente. a y b deben tener la misma cantidad de componentes."""
    a, b = _par(a, b)
    alfa, beta = numero(alfa, "α"), numero(beta, "β")
    return _limpia(alfa * a + beta * b, np.abs(alfa * a) + np.abs(beta * b))


def unitario(v):
    """v / ‖v‖: misma dirección y sentido, magnitud 1. El vector cero no tiene."""
    v = vector(v)
    if not np.any(v):
        raise ValueError("el vector cero no tiene dirección, así que no tiene vector unitario.")
    return v / _norma(v)


def desigualdad_triangulo(a, b, alfa=1, beta=1):
    """(‖α a + β b‖, |α| ‖a‖ + |β| ‖b‖, True si son iguales)."""
    c = combinacion(a, b, alfa, beta)
    a, b = _par(a, b)
    alfa, beta = numero(alfa, "α"), numero(beta, "β")
    izquierda = _norma(c)
    derecha = abs(alfa) * _norma(a) + abs(beta) * _norma(b)
    return izquierda, derecha, abs(derecha - izquierda) <= 1e-12 * derecha


def calculadora_suma(a_txt, b_txt, alfa_txt, beta_txt, n_cifras):
    try:
        a, b = _par(a_txt, b_txt)
        alfa, beta = numero(alfa_txt, "α"), numero(beta_txt, "β")
        c = combinacion(a, b, alfa, beta)
        izq, der, igual = desigualdad_triangulo(a, b, alfa, beta)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"a = {fmt(a, n_cifras)}    b = {fmt(b, n_cifras)}    α = {cifras(alfa, n_cifras)}    β = {cifras(beta, n_cifras)}")
    print(f"α a + β b = {fmt(c, n_cifras)}")
    print(f"‖α a + β b‖ = {cifras(izq, n_cifras)}")
    try:
        print(f"vector unitario = {fmt(unitario(c), n_cifras)}   (misma dirección y sentido, magnitud 1)")
    except ValueError as err:
        print("vector unitario:", err)
    print(f"‖α a + β b‖ = {cifras(izq, n_cifras)}  ≤  |α| ‖a‖ + |β| ‖b‖ = {cifras(der, n_cifras)}: "
          "se cumple la desigualdad del triángulo")
    print("   y con igualdad: α a y β b apuntan en el mismo sentido (o alguno es el vector cero)." if igual else
          "   con desigualdad estricta: α a y β b no apuntan en el mismo sentido.")
    if len(c) != 2:
        print("(El paralelogramo se dibuja solo para vectores del plano.)")
        return
    with np.errstate(all="ignore"):          # componentes cercanas a 1e100: matplotlib eleva al cuadrado
        _dibuja_paralelogramo(alfa * a, beta * b, c)


def _dibuja_paralelogramo(A, B, c):
    o = np.zeros(2)
    fig, ax = plt.subplots(figsize=(5.5, 4.5))
    for v, color, e in ((A, "tab:red", "α a"), (B, "tab:green", "β b"), (c, "navy", "α a + β b")):
        ax.quiver(*o, *v, angles="xy", scale_units="xy", scale=1, color=color, label=e, width=0.01)
    ax.plot([A[0], c[0]], [A[1], c[1]], "k:", lw=1)       # lados del paralelogramo
    ax.plot([B[0], c[0]], [B[1], c[1]], "k:", lw=1)
    _ejes_iguales(ax, [o, A, B, c])
    ax.grid(True); ax.legend(loc="best", fontsize=8)
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_title("Regla del paralelogramo: α a + β b")
    plt.show()


widgets.interact(calculadora_suma,
    a_txt=widgets.Text(value="300, 0", description="a =", continuous_update=False),
    b_txt=widgets.Text(value="0, 400", description="b =", continuous_update=False),
    alfa_txt=widgets.Text(value="1", description="α =", continuous_update=False),
    beta_txt=widgets.Text(value="1", description="β =", continuous_update=False),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras sig.", continuous_update=False));

# %%
# Casos de prueba de la sección 2 (resultado conocido)
PRUEBAS_2 = [
    ("(300, 0) + (0, 400) = (300, 400), magnitud 500 < 700; resta (300, -400); con α = -1, |α| ‖a‖ + ‖b‖ = 3",
     lambda: np.allclose(combinacion((300, 0), (0, 400), 1, 1), (300, 400))
             and cerca(desigualdad_triangulo((300, 0), (0, 400), 1, 1)[0], 500)
             and not desigualdad_triangulo((300, 0), (0, 400), 1, 1)[2]
             and np.allclose(combinacion((300, 0), (0, 400), 1, -1), (300, -400))
             and cerca(desigualdad_triangulo((1, 0), (2, 0), -1, 1)[1], 3)),
    ("unitario de (2, -1, 2) = (2/3, -1/3, 2/3)", lambda: np.allclose(unitario((2, -1, 2)), (2/3, -1/3, 2/3))),
    ("el vector (0, 0, 0) no tiene unitario", lambda: _rechaza(lambda: unitario((0, 0, 0)))),
    ("(3, 0) + (4, 0): igualdad 7 = 3 + 4",
     lambda: desigualdad_triangulo((3, 0), (4, 0), 1, 1)[2] and cerca(desigualdad_triangulo((3, 0), (4, 0), 1, 1)[0], 7)),
    ("un vector del plano con uno del espacio se rechaza", lambda: _rechaza(lambda: combinacion((1, 2), (1, 2, 3), 1, 1))),
]
for nombre, prueba in PRUEBAS_2:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 2).** Calcula a mano el vector unitario de $(1,-4,8)$: primero su magnitud, después divide cada componente entre ella. Escríbelo como texto, con fracciones o con decimales.

# %%
mi_unitario = None     # escribe tu vector como texto, por ejemplo: "1/3, 2/3, -2/3" o "0.3333, 0.6667, -0.6667"

if mi_unitario is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    u_ref = unitario((1, -4, 8))
    exacto = ", ".join(str(Fraction(c).limit_denominator(1000)) for c in u_ref)
    print(f"La calculadora da: ({exacto}) ≈", fmt(u_ref, 4))
    try:
        ok = bool(np.all(np.abs(vector(mi_unitario, 3, "tu vector") - u_ref) < 1e-3))
    except ValueError as err:
        print("Revisa tu valor:", err); ok = None
    if ok is not None:
        print("coinciden" if ok else
              "NO coinciden: la magnitud es la raíz de la suma de los cuadrados; divide cada componente entre ella")
