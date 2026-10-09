# ID: DIF-U4-NB04
# Notebook: dif/u4_derivada.ipynb · sección 4.4 notaciones y recta tangente
# Repositorio: dif/u4_derivada/04_recta_tangente.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 4. Notaciones $f'(x)$ y $dy/dx$; recta tangente
#
# Si $y=f(x)$, la derivada se escribe $f'(x)$, $y'$, $\dfrac{dy}{dx}$ o $\dfrac{d}{dx}f(x)$. Su valor en $a$ es $f'(a)$ o $\left.\dfrac{dy}{dx}\right|_{x=a}$. La recta tangente en $(a,f(a))$ es
# $$y=f(a)+f'(a)\,(x-a),$$
# y cerca de $a$ sirve para aproximar: $f(a+\Delta x)\approx f(a)+f'(a)\,\Delta x$.
#
# La calculadora da la ecuación de la tangente, compara la aproximación con el valor exacto en $a+\Delta x$ y dibuja las dos. Si los laterales del cociente tienden a $\pm\infty$ con el mismo signo, la tangente es vertical, $x=a$; si son distintos, no hay tangente.

# %%
def tangente(f, a):
    """(m, b) de la recta tangente y = m x + b en x = a; ValueError si no existe o es vertical."""
    izq, der, m = derivada_en(f, a)
    if m is None:
        infinitos = (sp.oo, -sp.oo)
        if izq is not None and der is not None and izq == der and izq in infinitos:
            raise ValueError(f"el cociente tiende a {a_texto(der)} por los dos lados: la tangente es vertical, x = {a}.")
        if (izq is None) != (der is None):
            lado, v = ("derecha", der) if izq is None else ("izquierda", izq)
            if v in infinitos:
                raise ValueError(f"x = {a} es un extremo del dominio y el cociente por la {lado} tiende a {a_texto(v)}: "
                                 f"la tangente es vertical, x = {a}, pero f'({a}) no existe.")
            raise ValueError(f"x = {a} es un extremo del dominio: solo existe la derivada por la {lado}, {a_texto(v)}.")
        raise ValueError(f"los laterales del cociente son {lado_texto(izq)} y {lado_texto(der)}: no hay recta tangente.")
    a = sp.nsimplify(a)
    return m, sp.simplify(f.subs(x, a) - m * a)


def calculadora_tangente(f_txt, a, dx, n_cifras):
    try:
        f = parsear(f_txt)
        m, b = tangente(f, a)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    a_, dx_ = sp.nsimplify(a), sp.nsimplify(dx)
    aprox, exacto = m * (a_ + dx_) + b, f.subs(x, a_ + dx_)
    print(f"f'({a}) = dy/dx en x = {a}: {m} ≈ {cifras(m, n_cifras)}")
    print(f"tangente: y = {cifras(m, n_cifras)} x + {cifras(b, n_cifras)}   (exacta: y = {sp.expand(m*x + b)})")
    print(f"en x = {a} + {dx}: tangente {cifras(aprox, n_cifras)} | f exacta {cifras(exacto, n_cifras)} "
          f"| error {cifras(sp.N(aprox - exacto), 3)}  ({n_cifras} cifras significativas)")
    fn = sp.lambdify(x, f, "numpy")
    xs = np.linspace(float(a) - 3 * max(abs(dx), 0.5), float(a) + 3 * max(abs(dx), 0.5), 400)
    with np.errstate(all="ignore"):
        ys = np.broadcast_to(fn(xs), xs.shape).astype(float)
    plt.figure(figsize=(5, 3.2)); plt.plot(xs, ys, color="gray", label="f")
    plt.plot(xs, float(m) * xs + float(b), "b--", label="tangente")
    plt.plot([a], [float(f.subs(x, a_))], "ro"); plt.grid(True); plt.legend(fontsize=8); plt.show()


widgets.interact(calculadora_tangente,
    f_txt=widgets.Text(value="1/x", description="f(x) ="),
    a=widgets.FloatText(value=2.0, description="a"), dx=widgets.FloatText(value=0.1, description="Δx"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 4 (resultado conocido)
PRUEBAS_4 = [
    ("x² en 3: y = 6x - 9", lambda: tangente(x**2, 3) == (6, -9)),
    ("1/x en 2: y = -x/4 + 1", lambda: tangente(1/x, 2) == (-sp.Rational(1, 4), 1)),
    ("√x en 4: y = x/4 + 1", lambda: tangente(sp.sqrt(x), 4) == (sp.Rational(1, 4), 1)),
    ("x³ en -1: y = 3x + 2", lambda: tangente(x**3, -1) == (3, 2)),
    ("Torricelli √(19.62 x) en 2: pendiente ≈ 1.5660",
     lambda: cerca(float(tangente(parsear("sqrt(19.62*x)"), 2)[0]), 1.566046, 1e-5)),
    ("|x| en 0: sin tangente (laterales -1 y 1)", lambda: _rechaza(lambda: tangente(sp.Abs(x), 0))),
    ("cbrt(x) en 0: tangente vertical, se avisa", lambda: _rechaza(lambda: tangente(parsear("cbrt(x)"), 0))),
    ("√x en 0: extremo con tangente vertical, se avisa", lambda: _rechaza(lambda: tangente(sp.sqrt(x), 0))),
]
for nombre, prueba in PRUEBAS_4:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 4).** Escribe a mano la recta tangente a $y=\sqrt{x}$ en $x=9$ en la forma $y=mx+b$. Compárala con la calculadora.

# %%
mi_m, mi_b = None, None      # escribe dos números, por ejemplo: 0.5, 1

ref = tangente(sp.sqrt(x), 9)
print("La calculadora da: m =", ref[0], ", b =", ref[1])
print("falta tu cálculo" if mi_m is None else
      ("coincide" if cerca(mi_m, float(ref[0])) and cerca(mi_b, float(ref[1])) else "NO coincide: y = f(9) + f'(9)(x - 9), con f'(9) = 1/(2·3)"))
