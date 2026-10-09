# ID: DIF-U6-NB02
# Notebook: dif/u6_aplicaciones.ipynb · sección 6.2 valor medio, diferenciales y Newton
# Repositorio: dif/u6_aplicaciones/02_valor_medio_newton.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Teorema del valor medio, diferenciales y método de Newton
#
# **Valor medio.** Si $f$ es continua en $[a,b]$ y derivable en $(a,b)$, existe $c$ en $(a,b)$ con $f'(c)=\dfrac{f(b)-f(a)}{b-a}$. **Diferencial:** $dy=f'(x)\,dx\approx\Delta y$.
#
# **Newton.** $x_{k+1}=x_k-\dfrac{f(x_k)}{f'(x_k)}$. Rápido cerca de una raíz simple; falla si $f'(x_k)=0$ y puede alejarse si se empieza lejos. La **bisección** de la Unidad 2 es más lenta pero siempre converge si $f$ es continua y cambia de signo en el intervalo. Newton se detiene cuando $|f(x_k)|$ y $|x_{k+1}-x_k|$ son menores que la tolerancia.

# %%
def newton(f, x0, tol=1e-10, max_iter=50):
    """(raíz, iteraciones, lista de x_k) por el método de Newton. ValueError si f'(x_k) = 0 o si no converge."""
    df = sp.diff(f, x).replace(sp.DiracDelta, lambda *args: sp.S.Zero)   # la derivada de sign(x) fuera de 0
    F, DF = sp.lambdify(x, f, "math"), sp.lambdify(x, df, "math")
    xk, historia = float(x0), [float(x0)]
    for k in range(1, max_iter + 1):
        d = DF(xk)
        if d == 0:
            raise ValueError(f"f'({cifras(xk, 6)}) = 0: la tangente es horizontal y no corta el eje. Cambia x0.")
        nuevo = xk - F(xk) / d
        historia.append(nuevo)
        if abs(nuevo) > 1e12:
            raise ValueError("los valores crecen sin cota: Newton no converge desde ese x0.")
        if abs(F(nuevo)) < tol and abs(nuevo - xk) < tol:
            return nuevo, k, historia
        xk = nuevo
    raise ValueError(f"no convergió en {max_iter} iteraciones: revisa x0 o usa bisección.")


def biseccion(f, a, b, tol=1e-10):
    """(raíz, iteraciones) por bisección; ValueError si f no cambia de signo en [a, b]."""
    F = sp.lambdify(x, f, "math")
    if F(a) * F(b) > 0:
        raise ValueError("f no cambia de signo en [a, b]: la bisección no garantiza una raíz ahí.")
    k = 0
    while b - a > tol:
        m = (a + b) / 2
        a, b = (a, m) if F(a) * F(m) <= 0 else (m, b)
        k += 1
    return (a + b) / 2, k


def valor_medio(f, a, b):
    """Puntos c en (a, b) con f'(c) = pendiente de la secante."""
    a, b = sp.nsimplify(a), sp.nsimplify(b)
    m = (f.subs(x, b) - f.subs(x, a)) / (b - a)
    sol = sp.solveset(sp.Eq(sp.diff(f, x), m), x, domain=sp.Interval.open(a, b))
    return m, sorted(sol, key=float)


def calculadora_newton(f_txt, x0, tol, max_iter, a, b, n_cifras):
    try:
        f = parsear(f_txt)
        if not 0 < tol < 1 or max_iter < 1:
            raise ValueError("la tolerancia debe estar entre 0 y 1 y el máximo de iteraciones ser al menos 1.")
        r, k, hist = newton(f, x0, tol, int(max_iter))
    except ValueError as err:
        print("Newton:", err)
    else:
        print("Newton:", " -> ".join(cifras(v, n_cifras) for v in hist[:8]), "..." if len(hist) > 8 else "")
        print(f"  raíz ≈ {cifras(r, n_cifras)} en {k} iteraciones ({n_cifras} cifras significativas)")
    try:
        rb, kb = biseccion(parsear(f_txt), a, b, tol)
        print(f"Bisección en [{a}, {b}]: raíz ≈ {cifras(rb, n_cifras)} en {kb} iteraciones")
    except ValueError as err:
        print("Bisección:", err)


widgets.interact(calculadora_newton,
    f_txt=widgets.Text(value="x**2 - 2", description="f(x) ="),
    x0=widgets.FloatText(value=1.0, description="x0"),
    tol=widgets.FloatText(value=1e-10, description="tolerancia"),
    max_iter=widgets.IntText(value=50, description="máx. iter."),
    a=widgets.FloatText(value=1.0, description="a (bisección)"), b=widgets.FloatText(value=2.0, description="b (bisección)"),
    n_cifras=widgets.IntSlider(value=8, min=1, max=12, description="cifras sig."));


def calculadora_valor_medio(f_txt, a, b, n_cifras):
    try:
        if not a < b:
            raise ValueError("el intervalo debe cumplir a < b.")
        m, cs = valor_medio(parsear(f_txt), a, b)
    except (ValueError, TypeError) as err:
        print("Revisa la entrada:", err); return
    print(f"pendiente de la secante en [{a}, {b}]: {cifras(m, n_cifras)}")
    print("puntos c con f'(c) igual a esa pendiente:", [cifras(c, n_cifras) for c in cs] or
          "ninguno: revisa si f es continua en [a, b] y derivable en (a, b)")


widgets.interact(calculadora_valor_medio,
    f_txt=widgets.Text(value="sqrt(x)", description="f(x) ="),
    a=widgets.FloatText(value=1.0, description="a"), b=widgets.FloatText(value=9.0, description="b"),
    n_cifras=widgets.IntSlider(value=6, min=1, max=10, description="cifras sig."));

# %%
# Casos de prueba de la sección 2 (resultado conocido)
PRUEBAS_2 = [
    ("Newton √2 desde 1 en 5 iteraciones", lambda: cerca(newton(x**2 - 2, 1)[0], math.sqrt(2), 1e-12)
                                                  and newton(x**2 - 2, 1)[1] == 5),
    ("bisección √2 en [1, 2]: 34 iteraciones", lambda: biseccion(x**2 - 2, 1, 2)[1] == 34),
    ("e^(-x) = x desde 1: 0.5671433", lambda: cerca(newton(sp.exp(-x) - x, 1)[0], 0.5671432904, 1e-9)),
    ("x³ - 2x - 5 desde 3: x1 = 2.36", lambda: cerca(newton(x**3 - 2*x - 5, 3)[2][1], 2.36, 1e-12)),
    ("x² - 4 desde 0: f'(0) = 0, se avisa", lambda: _rechaza(lambda: newton(x**2 - 4, 0))),
    ("cbrt(x) desde 1: no converge, se avisa", lambda: _rechaza(lambda: newton(parsear("cbrt(x)"), 1))),
    ("valor medio √x en [1, 9]: c = 4", lambda: valor_medio(sp.sqrt(x), 1, 9) == (sp.Rational(1, 4), [4])),
    ("valor medio x³ en [0, 2]: c = 2/√3", lambda: valor_medio(x**3, 0, 2)[1] == [2/sp.sqrt(3)]),
]
for nombre, prueba in PRUEBAS_2:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 2).** Haz a mano dos pasos de Newton para $x^3-2x-5=0$ desde $x_0=2$ y escribe $x_2$ con 5 cifras significativas.

# %%
mi_x2 = None             # escribe un número

ref = newton(x**3 - 2*x - 5, 2)[2][2]
if mi_x2 is None:                     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    print("La calculadora da x2 =", cifras(ref, 5))
    print("coincide" if abs(mi_x2 - ref) < 5e-5 else "NO coincide: x1 = 2 - f(2)/f'(2) = 2 - (-1)/10")
