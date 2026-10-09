# ID: DIF-U4-NB01
# Notebook: dif/u4_derivada.ipynb · sección 4.1 pendiente de una recta y de la secante
# Repositorio: dif/u4_derivada/01_pendiente_secante.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Pendiente de una recta y pendiente de la secante
#
# La pendiente de la recta que pasa por $(x_1,y_1)$ y $(x_2,y_2)$ es
# $$m=\frac{\Delta y}{\Delta x}=\frac{y_2-y_1}{x_2-x_1},\qquad x_1\neq x_2.$$
# En una recta da lo mismo con cualquier par de puntos. En una curva, la recta que pasa por $(a,f(a))$ y $(b,f(b))$ es una **secante** y su pendiente cambia con los puntos: es la **razón de cambio promedio** de $f$ entre $a$ y $b$.
#
# La calculadora recibe $f$, $a$ y $b$, da la pendiente exacta y redondeada, y dibuja la curva con la secante.

# %%
def calculadora_secante(f_txt, a, b, n_cifras):
    try:
        f = parsear(f_txt)
        m = pendiente_secante(f, a, b)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    fa, fb = f.subs(x, sp.nsimplify(a)), f.subs(x, sp.nsimplify(b))
    print(f"f(a) = {cifras(fa, n_cifras)}   f(b) = {cifras(fb, n_cifras)}")
    print(f"Δy/Δx = ({cifras(fb, n_cifras)} - {cifras(fa, n_cifras)}) / ({b} - {a})")
    print(f"pendiente de la secante = {m} ≈ {cifras(m, n_cifras)}  ({n_cifras} cifras significativas)")
    fn = sp.lambdify(x, f, "numpy")
    lo, hi = min(a, b), max(a, b); ancho = hi - lo
    xs = np.linspace(lo - ancho, hi + ancho, 400)
    with np.errstate(all="ignore"):
        ys = np.broadcast_to(fn(xs), xs.shape).astype(float)
    plt.figure(figsize=(5, 3.2)); plt.plot(xs, ys, color="gray", label="f")
    plt.plot(xs, float(fa) + float(m) * (xs - a), "b--", label=f"secante, m ≈ {cifras(m, 4)}")
    plt.plot([a, b], [float(fa), float(fb)], "ro"); plt.grid(True); plt.legend(fontsize=8); plt.show()


widgets.interact(calculadora_secante,
    f_txt=widgets.Text(value="x**2", description="f(x) ="),
    a=widgets.FloatText(value=1.0, description="a"), b=widgets.FloatText(value=3.0, description="b"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 1 (resultado conocido)
PRUEBAS_1 = [
    ("recta 0.5x + 1: pendiente 1/2 con (1, 3) y con (-2, 10)",
     lambda: pendiente_secante(parsear("0.5*x + 1"), 1, 3) == sp.Rational(1, 2)
             and pendiente_secante(parsear("0.5*x + 1"), -2, 10) == sp.Rational(1, 2)),
    ("x² de 1 a 3: pendiente 4", lambda: pendiente_secante(x**2, 1, 3) == 4),
    ("x² de 1 a 2: pendiente 3", lambda: pendiente_secante(x**2, 1, 2) == 3),
    ("costo 2000 + 15q + 0.05q² de 100 a 120: 26 pesos/pieza",
     lambda: pendiente_secante(parsear("2000 + 15*x + 0.05*x**2"), 100, 120) == 26),
    ("el orden de los puntos no cambia la pendiente",
     lambda: pendiente_secante(x**3, 1, 2) == pendiente_secante(x**3, 2, 1)),
    ("a = b se rechaza con mensaje", lambda: _rechaza(lambda: pendiente_secante(x**2, 2, 2))),
    ("punto fuera del dominio se rechaza", lambda: _rechaza(lambda: pendiente_secante(sp.sqrt(x), -1, 4))),
]
for nombre, prueba in PRUEBAS_1:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 1).** Calcula a mano la pendiente de la secante de $y=x^3$ entre $x=1$ y $x=2$. Escribe tu resultado y compáralo con la calculadora.

# %%
mi_pendiente = None       # escribe un número, por ejemplo: 2.5

ref = pendiente_secante(x**3, 1, 2)
print("La calculadora da:", ref)
print("falta tu cálculo" if mi_pendiente is None else
      ("coincide" if cerca(mi_pendiente, float(ref)) else "NO coincide: Δy = 2³ - 1³ y Δx = 2 - 1"))
