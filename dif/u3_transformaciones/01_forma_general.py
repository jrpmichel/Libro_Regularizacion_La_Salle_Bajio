# ID: DIF-U3-NB01
# Notebook: dif/u3_transformaciones.ipynb · sección 3.1 la forma general
# Repositorio: dif/u3_transformaciones/01_forma_general.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. La forma general $y=a\,f\big(b(x-h)\big)+k$
#
# Cada punto $(x_0,y_0)$ de la gráfica de $f$ tiene una imagen en la gráfica de $g$:
# $$(x_0,\,y_0)\ \longmapsto\ \left(h+\frac{x_0}{b},\ k+a\,y_0\right).$$
# Los parámetros de **afuera** ($a$, $k$) actúan sobre $y$ tal como se leen. Los de **adentro** ($b$, $h$) actúan sobre $x$: se divide entre $b$ y se suma $h$.
#
# La calculadora aplica la transformación a una función base del menú, muestra la imagen de tres puntos clave, el dominio y el rango, y dibuja las dos gráficas.

# %%
def calculadora_forma(base, a, b, h, k, n_cifras):
    try:
        f, puntos = BASES[base]
        g = transformar(f, a, b, h, k)
    except (KeyError, ValueError) as err:
        print("Revisa la entrada:", err); return
    print(f"f(x) = {f}      g(x) = {g}")
    print("Puntos clave (redondeados a", n_cifras, "cifras significativas):")
    for p in puntos:
        q = imagen(p, a, b, h, k)
        print(f"   ({cifras(p[0], n_cifras)}, {cifras(p[1], n_cifras)})  ->  ({cifras(q[0], n_cifras)}, {cifras(q[1], n_cifras)})"
              f"   comprobación g = {cifras(g.subs(x, q[0]), n_cifras)}")
    try:
        print("dominio de g:", conjunto_a_texto(continuous_domain(g, x, sp.S.Reals)))
        print("rango de g  :", conjunto_a_texto(function_range(g, x, continuous_domain(g, x, sp.S.Reals))))
    except Exception:
        print("sympy no pudo calcular el dominio o el rango en este caso.")
    fn, gn = sp.lambdify(x, f, "numpy"), sp.lambdify(x, g, "numpy")
    xs = np.linspace(-8, 8, 2001)
    with np.errstate(all="ignore"):
        yf = np.broadcast_to(fn(xs), xs.shape).astype(float)
        yg = np.broadcast_to(gn(xs), xs.shape).astype(float)
    yf = np.where(np.abs(yf) < 20, yf, np.nan); yg = np.where(np.abs(yg) < 20, yg, np.nan)
    plt.figure(figsize=(5, 3.2)); plt.plot(xs, yf, color="gray", label="f (base)")
    plt.plot(xs, yg, "b--", label="g (transformada)")
    plt.axhline(0, color="k", lw=0.6); plt.axvline(0, color="k", lw=0.6)
    plt.ylim(-10, 10); plt.grid(True); plt.legend(fontsize=8); plt.show()


widgets.interact(calculadora_forma,
    base=widgets.Dropdown(options=list(BASES), value="√x", description="base"),
    a=widgets.FloatSlider(value=-2, min=-4, max=4, step=0.5, description="a"),
    b=widgets.FloatSlider(value=2, min=-4, max=4, step=0.5, description="b"),
    h=widgets.FloatSlider(value=1, min=-5, max=5, step=0.5, description="h"),
    k=widgets.FloatSlider(value=3, min=-5, max=5, step=0.5, description="k"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 1 (resultado conocido)
def _img(base, a, b, h, k):
    return [imagen(p, a, b, h, k) for p in BASES[base][1]]

g31 = transformar(sp.sqrt(x), -2, 2, 1, 3)
PRUEBAS_1 = [
    ("√x con a=-2, b=2, h=1, k=3: (0,0)->(1,3), (1,1)->(3/2,1), (4,2)->(3,-1)",
     lambda: _img("√x", -2, 2, 1, 3) == [(1, 3), (sp.Rational(3, 2), 1), (3, -1)]),
    ("la imagen de cada punto está en la gráfica de g",
     lambda: all(sp.simplify(g31.subs(x, q[0]) - q[1]) == 0 for q in _img("√x", -2, 2, 1, 3))),
    ("dominio [1, ∞) y rango (-∞, 3]",
     lambda: continuous_domain(g31, x, sp.S.Reals) == sp.Interval(1, sp.oo)
             and function_range(g31, x, sp.Interval(1, sp.oo)) == sp.Interval(-sp.oo, 3)),
    ("x² con h=3: el vértice (0,0) va a (3,0)", lambda: _img("x²", 1, 1, 3, 0)[1] == (3, 0)),
    ("sen x con b=2: (π/2, 1) va a (π/4, 1)", lambda: _img("sen x", 1, 2, 0, 0)[1] == (sp.pi/4, 1)),
    ("1/x con h=3, k=2: dominio sin el 3",
     lambda: continuous_domain(transformar(1/x, 1, 1, 3, 2), x, sp.S.Reals)
             == sp.Union(sp.Interval.open(-sp.oo, 3), sp.Interval.open(3, sp.oo))),
    ("a = 0 se rechaza con mensaje", lambda: _rechaza(lambda: transformar(x**2, 0, 1, 0, 0))),
    ("b = 0 se rechaza con mensaje", lambda: _rechaza(lambda: transformar(x**2, 1, 0, 0, 0))),
]
for nombre, prueba in PRUEBAS_1:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 1).** Para $g(x)=2\,|x-1|+3$ (base $|x|$, con $a=2$, $b=1$, $h=1$, $k=3$), calcula a mano la imagen del punto $(-1,1)$ y comprueba que está en la gráfica de $g$. Escribe tu resultado y compáralo con la calculadora.

# %%
mi_imagen = None       # escribe una pareja (x, y), por ejemplo: (1, 2)

ref = imagen((-1, 1), 2, 1, 1, 3)
print("La calculadora da:", ref, "| g en ese x:", transformar(sp.Abs(x), 2, 1, 1, 3).subs(x, ref[0]))
print("falta tu cálculo" if mi_imagen is None else
      ("coincide" if all(cerca(float(u), float(v)) for u, v in zip(mi_imagen, ref)) else "NO coincide: ¿dividiste entre b y sumaste h?"))
