# ID: INT-U3-NB03
# Notebook: int/u3_antiderivada.ipynb · sección 3.3 condiciones iniciales
# Repositorio: int/u3_antiderivada/03_condiciones.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Condiciones iniciales
#
# * Primer orden: $y'=f(x)$, $y(x_0)=y_0$. Solución $y=F(x)+C$ con $C=y_0-F(x_0)$.
# * Segundo orden: $y''=g(x)$, $y'(x_0)=v_0$, $y(x_0)=y_0$. Se integra dos veces y cada constante se fija **en cuanto aparece**.
#
# La calculadora muestra cada constante despejada, la solución y su gráfica. Rechaza un $x_0$ donde $f$ no está definida y recuerda que la solución vale en el intervalo que contiene a $x_0$.

# %%
def resolver_pvi(f, x0, y0, orden=1, v0=0):
    """(y, y', constantes) del problema de valor inicial de primer o segundo orden."""
    x0 = sp.nsimplify(x0); y0 = sp.nsimplify(y0); v0 = sp.nsimplify(v0)
    valor_real(f, x0)                                   # ValueError si f no existe en x0
    if orden == 1:
        F = antiderivada(f)
        C = sp.simplify(y0 - F.subs(x, x0))
        return _ordenar(F + C), f, {"C": C}
    G = antiderivada(f)
    C1 = sp.simplify(v0 - G.subs(x, x0))
    yp = G + C1
    Y = antiderivada(yp) if yp.has(x) else yp * x
    C2 = sp.simplify(y0 - Y.subs(x, x0))
    return _ordenar(Y + C2), _ordenar(yp), {"C1": C1, "C2": C2}


def _ordenar(e):
    """Polinomios desarrollados; lo demás, simplificado."""
    return sp.expand(e) if e.is_polynomial(x) else sp.simplify(e)


def calculadora_pvi(f_txt, orden, x0, y0, v0, x_eval, n_cifras):
    try:
        f = parsear(f_txt)
        y, yp, cons = resolver_pvi(f, x0, y0, 1 if orden.startswith("1") else 2, v0)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    for nombre, valor in cons.items():
        print(f"{nombre} = {valor}")
    print("y(x) =", y)
    if not orden.startswith("1"):
        print("y'(x) =", yp)
    try:
        print(f"y({x_eval}) ≈ {cifras(valor_real(y, x_eval), n_cifras)}   ({n_cifras} cifras significativas, redondeado)")
    except ValueError as err:
        print("No se puede evaluar en ese punto:", err)
    print("La solución vale en el intervalo que contiene a x0 donde f es continua.")
    xs = np.linspace(float(x0) - 3, float(x0) + 3, 400)
    with np.errstate(all="ignore"):
        fig, ax = plt.subplots(figsize=(5.5, 3))
        ax.plot(xs, numerica(y)(xs), color="navy", lw=2, label="y(x)")
        if not orden.startswith("1"):
            ax.plot(xs, numerica(yp)(xs), "--", color="red", label="y'(x)")
    ax.plot([float(x0)], [float(y0)], "o", color="red")
    ax.axhline(0, color="black", lw=0.6); ax.set_xlabel("x"); ax.legend(); plt.show()


widgets.interact(calculadora_pvi,
    f_txt=widgets.Text(value="3*x**2 - 3", description="f o g ="),
    orden=widgets.Dropdown(options=["1: y' = f", "2: y'' = g"], value="1: y' = f", description="orden"),
    x0=widgets.FloatText(value=2, description="x0"), y0=widgets.FloatText(value=5, description="y0"),
    v0=widgets.FloatText(value=0, description="y'(x0)"),
    x_eval=widgets.FloatText(value=3, description="evaluar en"),
    n_cifras=widgets.IntSlider(value=5, min=1, max=10, description="cifras sig."));

# %%
# Casos de prueba de la sección 3 (resultado conocido)
PRUEBAS_3 = [
    ("y' = 2x, y(0) = 1 -> x^2 + 1", lambda: misma_funcion(resolver_pvi(2 * x, 0, 1)[0], x**2 + 1)),
    ("y'' = -9.8, y'(0) = 12, y(0) = 1.5 -> y(1) = 8.6",
     lambda: cerca(valor_real(resolver_pvi(sp.Rational(-49, 5) + 0 * x, 0, 1.5, 2, 12)[0], 1), 8.6)),
    ("y' = 1/x, y(1) = 2 -> y(e) = 3", lambda: cerca(valor_real(resolver_pvi(1 / x, 1, 2)[0], math.e), 3)),
    ("y' = 1/x con x0 = 0 se rechaza", lambda: _rechaza(lambda: resolver_pvi(1 / x, 0, 1))),
    ("y' = cos x, y(0) = 0 -> sen x", lambda: misma_funcion(resolver_pvi(sp.cos(x), 0, 0)[0], sp.sin(x))),
    ("y'' = 6x, y'(0) = 0, y(0) = 0 -> x^3", lambda: misma_funcion(resolver_pvi(6 * x, 0, 0, 2, 0)[0], x**3)),
]
for nombre, prueba in PRUEBAS_3:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 3).** Resuelve a mano $y'=3x^2$ con $y(1)=5$.

# %%
mi_y = None     # escribe tu solución como texto, por ejemplo: "x**3 + 1"

if mi_y is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    ref = resolver_pvi(3 * x**2, 1, 5)[0]
    print("La calculadora da: y =", ref)
    try:
        ok = misma_funcion(parsear(str(mi_y)), ref)
    except ValueError as err:
        print(err); ok = False
    print("coinciden" if ok else "NO coinciden: después de integrar, sustituye x = 1 y despeja C para que y(1) = 5")
