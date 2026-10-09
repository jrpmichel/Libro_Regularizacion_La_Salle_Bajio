# ID: DIF-U2-NB03
# Notebook: dif/u2_limite.ipynb · sección 2.3 límites en el infinito
# Repositorio: dif/u2_limite/03_infinito.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Límites en el infinito
#
# $\lim_{x\to\infty}f(x)=L$ significa que $f(x)$ se acerca a $L$ cuando $x$ crece sin tope; entonces la recta $y=L$ es una **asíntota horizontal**. En $-\infty$ se mira hacia la izquierda, y el resultado puede ser otro: cuidado con las raíces, porque $\sqrt{x^2}=|x|=-x$ cuando $x<0$.
#
# La calculadora tabula $f$ en $x=\pm10,\ \pm100,\dots$, calcula el límite con `sympy` y dibuja la función con su asíntota.

# %%
def limite_infinito(texto, direccion="+∞"):
    """Límite de f cuando x -> +∞ o x -> -∞."""
    if direccion not in ("+∞", "-∞"):
        raise ValueError("La dirección debe ser '+∞' o '-∞'.")
    f = parsear(texto)
    return sp.limit(f, x, sp.oo if direccion == "+∞" else -sp.oo)


def calculadora_infinito(texto, direccion, n_cifras):
    try:
        L = limite_infinito(texto, direccion)
        f = sp.lambdify(x, parsear(texto), "numpy")
    except ValueError as err:
        print("Revisa la entrada:", err); return
    s = 1 if direccion == "+∞" else -1
    print(f"f(x) = {texto}   cuando x -> {direccion}   (valores redondeados a {n_cifras} cifras significativas)")
    for k in range(1, 7):
        v = s * 10.0**k
        print(f"   x = {v:>10.0e}   f(x) = {cifras(evaluar(f, v), n_cifras)}")
    print("límite:", a_texto(L))
    if L.is_finite:
        print(f"asíntota horizontal: y = {a_texto(L)} (la curva puede cruzarla en algún x finito).")
    else:
        print("no hay asíntota horizontal de ese lado.")
    xs = np.linspace(-50, 50, 2001)
    with np.errstate(all="ignore"):
        ys = np.broadcast_to(f(xs), xs.shape).astype(float)
    ys = np.where(np.abs(ys) < 1e3, ys, np.nan)
    plt.figure(figsize=(5, 3)); plt.plot(xs, ys, "b-")
    if L.is_finite:
        plt.axhline(float(L), color="gray", ls="--", label=f"y = {a_texto(L)}"); plt.legend(fontsize=8)
    plt.grid(True); plt.xlabel("x"); plt.ylabel("f(x)"); plt.show()


widgets.interact(calculadora_infinito,
    texto=widgets.Text(value="(3*x**2 + 2*x)/(x**2 + 4)", description="f(x) ="),
    direccion=widgets.Dropdown(options=["+∞", "-∞"], value="+∞", description="x ->"),
    n_cifras=widgets.IntSlider(value=6, min=1, max=10, description="cifras sig."));

# %% [markdown]
# **Valor final de una carga.** Un capacitor que se carga a través de una resistencia sigue $v(t)=V_f\,(1-e^{-t/\tau})$. Mueve $\tau$ y $V_f$ y mira cuánto falta para el valor final en $t=\tau,\ 3\tau,\ 5\tau$: la fracción no depende de $\tau$.

# %%
def carga(Vf, tau):
    if tau <= 0:
        print("Revisa la entrada: la constante de tiempo debe ser positiva."); return
    for k in (1, 3, 5):
        v = Vf * (1 - math.exp(-k))
        print(f"t = {k}τ = {cifras(k*tau, 4)} s:  v = {cifras(v, 4)} V  ({cifras(100*v/Vf, 4)} % del valor final)")
    t = np.linspace(0, 6*tau, 300)
    plt.figure(figsize=(5, 3)); plt.plot(t, Vf*(1 - np.exp(-t/tau)), "b-"); plt.axhline(Vf, color="gray", ls="--")
    plt.grid(True); plt.xlabel("t (s)"); plt.ylabel("v (V)"); plt.show()


widgets.interact(carga, Vf=widgets.FloatText(value=12.0, description="V_f (V)"),
                 tau=widgets.FloatText(value=0.5, description="τ (s)"));

# %%
# Casos de prueba de la sección 3 (resultado conocido)
PRUEBAS_3 = [
    ("(3x²+2x)/(x²+4) -> 3 en +∞ y en -∞",
     lambda: limite_infinito("(3*x**2 + 2*x)/(x**2 + 4)") == 3 and limite_infinito("(3*x**2 + 2*x)/(x**2 + 4)", "-∞") == 3),
    ("(5x+1)/(x²-2) -> 0", lambda: limite_infinito("(5*x + 1)/(x**2 - 2)") == 0),
    ("x³/(x²+1) -> +∞ (sin asíntota horizontal)", lambda: limite_infinito("x**3/(x**2 + 1)") == sp.oo),
    ("(2x+1)/sqrt(x²+1) -> 2 en +∞ y -2 en -∞",
     lambda: limite_infinito("(2*x + 1)/sqrt(x**2 + 1)") == 2 and limite_infinito("(2*x + 1)/sqrt(x**2 + 1)", "-∞") == -2),
    ("12(1 - e^(-x/0.5)) -> 12", lambda: limite_infinito("12*(1 - exp(-x/0.5))") == 12),
    ("sqrt(x²+3x) - x -> 3/2 (no 0)", lambda: limite_infinito("sqrt(x**2 + 3*x) - x") == sp.Rational(3, 2)),
    ("una dirección inválida se rechaza", lambda: _rechaza(lambda: limite_infinito("x", "infinito"))),
]
for nombre, prueba in PRUEBAS_3:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 3).** Calcula a mano $\displaystyle\lim_{x\to\infty}\frac{4x^2-1}{2x^2+x}$ dividiendo numerador y denominador entre $x^2$. Escribe tu resultado y compáralo con la calculadora. Anota qué término domina en cada parte de la fracción.

# %%
mi_limite = None        # por ejemplo: 2

L = limite_infinito("(4*x**2 - 1)/(2*x**2 + x)")
print("La calculadora da:", L)
print("falta tu cálculo" if mi_limite is None else ("coincide" if cerca(mi_limite, float(L)) else "NO coincide: revisa qué potencia de x domina"))
