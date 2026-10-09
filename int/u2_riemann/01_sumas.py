# ID: INT-U2-NB01
# Notebook: int/u2_riemann.ipynb · sección 2.1 sumas izquierda, derecha y de punto medio
# Repositorio: int/u2_riemann/01_sumas.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Suma izquierda, derecha y de punto medio
#
# Con $\Delta x=\frac{b-a}{n}$ y $x_k=a+k\,\Delta x$:
#
# $$L_n=\sum_{k=0}^{n-1}f(x_k)\,\Delta x,\qquad R_n=\sum_{k=1}^{n}f(x_k)\,\Delta x,\qquad M_n=\sum_{k=1}^{n}f\!\left(\tfrac{x_{k-1}+x_k}{2}\right)\Delta x.$$
#
# | Si en $[a,b]$… | entonces |
# |---|---|
# | $f$ crece | $L_n\le A\le R_n$ |
# | $f$ decrece | $R_n\le A\le L_n$ |
# | $f''>0$ (cóncava hacia arriba) | $M_n\le A$ |
# | $f''<0$ (cóncava hacia abajo) | $M_n\ge A$ |
#
# La calculadora revisa el signo de $f'$ y $f''$ en 400 puntos del intervalo y dice qué suma sobreestima. Si el signo cambia, no concluye.

# %%
def sumas_lrm(f, a, b, m):
    """(L_m, R_m, M_m) de la expresión f en [a, b] con m franjas iguales."""
    if not a < b:
        raise ValueError("el intervalo debe cumplir a < b.")
    if int(m) != m or m < 1:
        raise ValueError("el número de franjas debe ser un entero positivo.")
    m = int(m)
    for c in np.linspace(a, b, min(m, 50) * 2 + 1):          # detecta puntos fuera del dominio
        valor_real(f, c)
    F = numerica(f)
    xs = np.linspace(a, b, m + 1)
    dx = (b - a) / m
    return float(F(xs[:-1]).sum() * dx), float(F(xs[1:]).sum() * dx), float(F((xs[:-1] + xs[1:]) / 2).sum() * dx)


def _signo_constante(expr, a, b):
    """+1 o -1 si expr no cambia de signo en (a, b); 0 si es idénticamente cero; None si cambia de signo."""
    v = numerica(expr)(np.linspace(a, b, 402)[1:-1])
    if np.allclose(v, 0, atol=1e-12): return 0
    if np.all(v > 0): return 1
    if np.all(v < 0): return -1
    if np.all(v >= 0): return 1      # se anula en puntos aislados (como x^2 en 0)
    if np.all(v <= 0): return -1
    return None


def clasificar_sumas(f, a, b):
    """Texto con qué suma sobreestima y cuál subestima, según monotonía y concavidad."""
    s1, s2 = _signo_constante(sp.diff(f, x), a, b), _signo_constante(sp.diff(f, x, 2), a, b)
    partes = []
    if s1 == 0: partes.append("f es constante: L_n, R_n y M_n son exactas")
    elif s1 == 1: partes.append("f crece: L_n subestima y R_n sobreestima")
    elif s1 == -1: partes.append("f decrece: L_n sobreestima y R_n subestima")
    else: partes.append("f no es monótona en [a, b]: no se puede ordenar L_n y R_n sin más información")
    if s2 == 0 and s1 != 0: partes.append("f'' = 0 (f es una recta): M_n es exacta")
    elif s2 == 1: partes.append("f'' > 0: M_n subestima")
    elif s2 == -1: partes.append("f'' < 0: M_n sobreestima")
    elif s2 is None: partes.append("la concavidad cambia: no se puede decir si M_n sobreestima o subestima")
    return "; ".join(partes)


def calculadora_sumas(f_txt, a, b, m, regla, n_cifras):
    try:
        f = parsear(f_txt)
        L, R, M = sumas_lrm(f, a, b, m)
    except (ValueError, TypeError) as err:
        print("Revisa la entrada:", err); return
    print(f"L_{m} ≈ {cifras(L, n_cifras)}   R_{m} ≈ {cifras(R, n_cifras)}   M_{m} ≈ {cifras(M, n_cifras)}"
          f"   ({n_cifras} cifras significativas, redondeado)")
    print(clasificar_sumas(f, a, b))
    F = numerica(f)
    xs = np.linspace(a, b, m + 1); dx = (b - a) / m
    pos = {"izquierda": xs[:-1], "derecha": xs[1:], "punto medio": (xs[:-1] + xs[1:]) / 2}[regla]
    fig, ax = plt.subplots(figsize=(5.5, 3))
    for x0, xm in zip(xs[:-1], pos):
        ax.add_patch(plt.Rectangle((x0, 0), dx, F(np.array([xm]))[0], facecolor="0.85", edgecolor="black"))
    ax.plot(pos, F(pos), "o", color="red")
    xx = np.linspace(a, b, 300); ax.plot(xx, F(xx), color="navy"); ax.axhline(0, color="black", lw=0.7)
    ax.set_xlabel("x"); ax.set_ylabel("f(x)")
    ax.set_title(f"suma {regla}, n = {m}"); plt.show()


widgets.interact(calculadora_sumas,
    f_txt=widgets.Text(value="x**2", description="f(x) ="),
    a=widgets.FloatText(value=0, description="a"), b=widgets.FloatText(value=1, description="b"),
    m=widgets.IntSlider(value=4, min=1, max=60, description="n"),
    regla=widgets.Dropdown(options=["izquierda", "derecha", "punto medio"], value="punto medio", description="regla"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 1 (resultado conocido)
_L, _R, _M = sumas_lrm(x**2, 0, 1, 4)
PRUEBAS_1 = [
    ("x^2 en [0, 1], n = 4: L = 0.21875, R = 0.46875, M = 0.328125",
     lambda: cerca(_L, 0.21875) and cerca(_R, 0.46875) and cerca(_M, 0.328125)),
    ("1/x en [1, 2], n = 4: M ≈ 0.691220", lambda: cerca(sumas_lrm(1 / x, 1, 2, 4)[2], 0.6912198912, 1e-9)),
    ("constante 2 en [0, 5]: las tres dan 10", lambda: all(cerca(v, 10) for v in sumas_lrm(2 + 0 * x, 0, 5, 3))),
    ("recta 5x + 2 en [0, 2], n = 3: M exacto = 14", lambda: cerca(sumas_lrm(5 * x + 2, 0, 2, 3)[2], 14)),
    ("hidrograma 30/(1 + x): decrece y f'' > 0", lambda: clasificar_sumas(30 / (1 + x), 0, 4).startswith("f decrece")
                                                     and "M_n subestima" in clasificar_sumas(30 / (1 + x), 0, 4)),
    ("recta 4x - 3: f'' = 0 y M_n exacta", lambda: "M_n es exacta" in clasificar_sumas(4 * x - 3, 1, 2)),
    ("a >= b, n = 0 o f fuera del dominio se rechazan", lambda: _rechaza(lambda: sumas_lrm(x, 1, 0, 4))
                                                              and _rechaza(lambda: sumas_lrm(x, 0, 1, 0))
                                                              and _rechaza(lambda: sumas_lrm(sp.sqrt(x), -1, 1, 4))),
]
for nombre, prueba in PRUEBAS_1:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 1).** Calcula a mano $L_3$, $R_3$ y $M_3$ para $f(x)=2x+1$ en $[0,3]$.

# %%
mi_L3 = None    # escribe un número
mi_R3 = None    # escribe un número
mi_M3 = None    # escribe un número

if None in (mi_L3, mi_R3, mi_M3):       # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    ref = sumas_lrm(2 * x + 1, 0, 3, 3)
    print("La calculadora da: L_3 =", cifras(ref[0], 4), " R_3 =", cifras(ref[1], 4), " M_3 =", cifras(ref[2], 4))
    print("coinciden" if all(cerca(u, v, 1e-6) for u, v in zip((mi_L3, mi_R3, mi_M3), ref)) else
          "NO coinciden: con Δx = 1, las alturas izquierdas son f(0), f(1), f(2) y las de punto medio f(0.5), f(1.5), f(2.5)")
