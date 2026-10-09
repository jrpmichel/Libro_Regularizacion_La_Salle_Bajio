# ID: DIF-U2-NB06
# Notebook: dif/u2_limite.ipynb · sección 2.6 continuidad y bisección
# Repositorio: dif/u2_limite/06_continuidad.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6. Continuidad, tipos de discontinuidad y método de bisección
#
# $f$ es continua en $a$ si se cumplen tres condiciones: $f(a)$ está definida, $\lim_{x\to a}f(x)$ existe y los dos coinciden. Si falla alguna, la discontinuidad se clasifica con los límites laterales:
#
# - **evitable**: los laterales coinciden y son finitos, pero $f(a)$ no existe o es otro número;
# - **de salto**: los laterales son finitos y distintos;
# - **infinita**: algún lateral es $\pm\infty$.
#
# La primera calculadora clasifica el punto frontera de una función por partes. La segunda aplica el **método de bisección**, que se apoya en el teorema del valor intermedio y por eso exige que $f$ sea continua en todo $[a,b]$. Antes de iterar, la calculadora lo comprueba.

# %%
def clasificar(izq_txt, der_txt, a, valor_txt):
    """'continua', 'evitable', 'de salto' o 'infinita', con los laterales."""
    li, ld = laterales(izq_txt, der_txt, a)
    valor = None if str(valor_txt).strip().lower() in ("no", "") else numero(valor_txt)
    if li.is_finite and ld.is_finite:
        if sp.simplify(li - ld) != 0:
            return "de salto", li, ld
        return ("continua" if valor is not None and sp.simplify(valor - li) == 0 else "evitable"), li, ld
    return "infinita", li, ld


def calculadora_continuidad(izq_txt, der_txt, a, valor_txt):
    try:
        tipo, li, ld = clasificar(izq_txt, der_txt, a, valor_txt)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"x < {a}: {izq_txt}    x > {a}: {der_txt}    f({a}) = {valor_txt}")
    print(f"izquierda {a_texto(li)} | derecha {a_texto(ld)}  ->  {tipo.upper()}")
    if tipo == "de salto":
        print("tamaño del salto (derecha - izquierda):", a_texto(ld - li))


widgets.interact(calculadora_continuidad,
    izq_txt=widgets.Text(value="0", description="izquierda"),
    der_txt=widgets.Text(value="10 + 18*x", description="derecha"),
    a=widgets.Text(value="0", description="a ="),
    valor_txt=widgets.Text(value="0", description="f(a) ="));

# %%
def biseccion(texto, a, b, tol=1e-6, max_pasos=100):
    """(raíz aproximada, pasos, filas) o ValueError si no se cumplen las hipótesis."""
    f_expr = parsear(texto)
    a, b = float(numero(a)), float(numero(b))
    if not a < b:
        raise ValueError("Necesitas a < b.")
    if not tol > 0:
        raise ValueError("La tolerancia debe ser positiva.")
    dominio = continuous_domain(f_expr, x, sp.Interval(sp.nsimplify(a), sp.nsimplify(b)))
    if dominio != sp.Interval(sp.nsimplify(a), sp.nsimplify(b)):
        raise ValueError(f"f no es continua en todo [{cifras(a)}, {cifras(b)}] (es continua en {conjunto_a_texto(dominio)}). "
                         "El teorema del valor intermedio no aplica: un cambio de signo puede venir de un polo.")
    f = sp.lambdify(x, f_expr, "math")
    fa, fb = f(a), f(b)
    if fa * fb > 0:
        raise ValueError(f"f(a) = {cifras(fa)} y f(b) = {cifras(fb)} tienen el mismo signo: la bisección no puede empezar.")
    filas, pasos = [], 0
    while b - a > tol and pasos < max_pasos:
        c = (a + b) / 2
        filas.append((pasos, a, b, c, f(c)))
        if f(a) * f(c) <= 0: b = c
        else: a = c
        pasos += 1
    return (a + b) / 2, pasos, filas


def calculadora_biseccion(texto, a, b, tol_exp, n_cifras):
    try:
        c, n, filas = biseccion(texto, a, b, 10.0**tol_exp)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    f = sp.lambdify(x, parsear(texto), "math")
    print(f"{'paso':>4} | {'a':>12} | {'b':>12} | {'c':>12} | {'f(c)':>12}")
    for p, a_, b_, c_, fc in filas[:8]:
        print(f"{p:>4} | {cifras(a_, n_cifras):>12} | {cifras(b_, n_cifras):>12} | {cifras(c_, n_cifras):>12} | {cifras(fc, n_cifras):>12}")
    if len(filas) > 8:
        print("   ...")
    print(f"raíz ≈ {cifras(c, n_cifras)} tras {n} pasos;  f(raíz) = {f(c):.1e} (debe ser pequeño)")


widgets.interact(calculadora_biseccion,
    texto=widgets.Text(value="x**3 - x - 2", description="f(x) ="),
    a=widgets.Text(value="1", description="a ="), b=widgets.Text(value="2", description="b ="),
    tol_exp=widgets.IntSlider(value=-6, min=-12, max=-1, description="tol = 10^"),
    n_cifras=widgets.IntSlider(value=7, min=1, max=12, description="cifras sig."));

# %%
# Casos de prueba de la sección 6 (resultado conocido)
PRUEBAS_6 = [
    ("x³ - x - 2 en [1, 2]: raíz 1.5214 en 20 pasos",
     lambda: (lambda r: round(r[0], 4) == 1.5214 and r[1] == 20)(biseccion("x**3 - x - 2", "1", "2"))),
    ("cos x - x en [0, 1]: raíz 0.7391 (radianes)", lambda: round(biseccion("cos(x) - x", "0", "1")[0], 4) == 0.7391),
    ("x/(x-2) - 3 en [0, 2.5] se rechaza: polo en x = 2 (Laboratorio 2.4)",
     lambda: _rechaza(lambda: biseccion("x/(x - 2) - 3", "0", "2.5"))),
    ("x/(x-2) - 3 en [2.5, 4]: raíz 3", lambda: round(biseccion("x/(x - 2) - 3", "2.5", "4")[0], 5) == 3.0),
    ("x² + 1 en [-1, 1]: sin cambio de signo, se rechaza", lambda: _rechaza(lambda: biseccion("x**2 + 1", "-1", "1"))),
    ("controlador P(e) en e = 0: salto de 10", lambda: (lambda r: r[0] == "de salto" and r[2] - r[1] == 10)(clasificar("0", "10 + 18*x", "0", "0"))),
    ("controlador P(e) en e = 5: continua", lambda: clasificar("10 + 18*x", "100", "5", "100")[0] == "continua"),
    ("(x²-1)/(x²-3x+2): evitable en 1 e infinita en 2",
     lambda: clasificar("(x**2 - 1)/(x**2 - 3*x + 2)", "(x**2 - 1)/(x**2 - 3*x + 2)", "1", "no")[0] == "evitable"
             and clasificar("(x**2 - 1)/(x**2 - 3*x + 2)", "(x**2 - 1)/(x**2 - 3*x + 2)", "2", "no")[0] == "infinita"),
]
for nombre, prueba in PRUEBAS_6:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 6).** Haz a mano dos pasos de bisección para $f(x)=x^3-x-2$ en $[1,2]$: evalúa $f$ en el punto medio, decide qué mitad conservas y repite. Escribe el intervalo que te queda y compáralo con las dos primeras filas de la calculadora.

# %%
mi_intervalo = None     # por ejemplo: (1.5, 1.75)

_, _, filas = biseccion("x**3 - x - 2", "1", "2")
ref = (filas[2][1], filas[2][2])
print("La calculadora, después de dos pasos, conserva:", ref)
print("falta tu cálculo" if mi_intervalo is None else
      ("coincide" if all(cerca(u, v) for u, v in zip(mi_intervalo, ref)) else "NO coincide: revisa el signo de f en el punto medio"))
