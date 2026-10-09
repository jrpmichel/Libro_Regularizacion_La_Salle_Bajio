# ID: INT-U3-NB01
# Notebook: int/u3_antiderivada.ipynb · sección 3.1 integral indefinida y constante de integración
# Repositorio: int/u3_antiderivada/01_indefinida.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Integral indefinida y constante de integración
#
# $F$ es antiderivada de $f$ en un intervalo si $F'(x)=f(x)$. Todas las antiderivadas difieren en una constante:
#
# $$\int f(x)\,dx=F(x)+C,\qquad \int x^n\,dx=\frac{x^{n+1}}{n+1}+C\quad(n\neq-1).$$
#
# La calculadora pide a `sympy` una antiderivada, le pone el valor absoluto al logaritmo cuando hace falta ($\ln|x|$), **siempre** escribe la constante $C$ y comprueba que $F'=f$. El control deslizante mueve $C$: la curva sube o baja, pero su pendiente en cada $x$ no cambia.

# %%
_PERMITIDAS = (sp.exp, sp.log, sp.sin, sp.cos, sp.tan, sp.sec, sp.Abs, sp.atan, sp.asin, sp.acos)
_PUNTOS = [-2.9, -2.1, -1.3, -0.7, 0.7, 1.3, 2.1, 2.9]


def con_valor_absoluto(F):
    """ln(u) -> ln|u| salvo que u sea siempre positivo: la antiderivada de u'/u vale donde u > 0 y donde u < 0."""
    return F.replace(sp.log, lambda a: sp.log(a) if (a.is_positive or a.has(sp.Abs)) else sp.log(sp.Abs(a)))


def sin_abs_en_log(f):
    """ln|u| -> ln u, solo para integrar (la derivada de las dos es u'/u); el valor absoluto se repone al final."""
    return f.replace(lambda e: isinstance(e, sp.log) and e.args[0].has(sp.Abs),
                     lambda e: sp.log(e.args[0].replace(sp.Abs, lambda a: a)))


def misma_funcion(g, h):
    """True si g y h coinciden: primero con sympy; si no decide, en los puntos donde las dos son reales."""
    d = sp.simplify(g - h)
    if d == 0:
        return True
    comparados = 0
    for c in _PUNTOS:
        try:
            if abs(valor_real(g, c) - valor_real(h, c)) > 1e-9:
                return False
            comparados += 1
        except (ValueError, TypeError, ZeroDivisionError):
            continue
    return comparados >= 3


def antiderivada(f):
    """Una antiderivada de f (con C = 0), simplificada y con ln|..| donde corresponde."""
    F = sp.integrate(sin_abs_en_log(f), x)
    if F.has(sp.Integral):
        raise ValueError("sympy no encontró una antiderivada para esta función.")
    raras = {type(a).__name__ for a in F.atoms(sp.Function) if not isinstance(a, _PERMITIDAS)}
    if raras:
        raise ValueError(f"la antiderivada usa funciones especiales ({', '.join(sorted(raras))}): no es elemental. "
                         "En la Unidad 4 verás cómo integrar estos casos con métodos numéricos.")
    F = sp.expand(F) if F.is_polynomial(x) else con_valor_absoluto(sp.trigsimp(sp.simplify(F)))
    if not misma_funcion(derivada_legible(F), sin_abs_en_log(f)):
        raise ValueError("la comprobación F' = f falló; revisa la función.")
    return F


def derivada_legible(F):
    """F' sin el Piecewise que sympy agrega al derivar ln|u| (la derivada de ln|u| es u'/u)."""
    G = F.replace(lambda e: isinstance(e, sp.log) and e.args[0].has(sp.Abs),
                  lambda e: sp.log(e.args[0].replace(sp.Abs, lambda a: a)))
    return sp.simplify(sp.diff(G, x))


def aviso_potencia(f):
    """Mensaje si f tiene un término x**(-1), donde la regla de la potencia no aplica."""
    terminos = sp.Add.make_args(sp.expand(f))
    if any(sp.simplify(t * x).is_number for t in terminos if t.has(x)):
        return "Aviso: la regla de la potencia no aplica a x^(-1) (pediría dividir entre 0); la tabla da ln|x|."
    return ""


def calculadora_indefinida(f_txt, C):
    try:
        f = parsear(f_txt)
        F = antiderivada(f)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"∫ ({f}) dx = {F} + C")
    print("Comprobación: d/dx [", F, "] =", derivada_legible(F))
    if aviso_potencia(f):
        print(aviso_potencia(f))
    Fn = numerica(F)
    xs = np.linspace(-3, 3, 601)
    with np.errstate(all="ignore"):
        fig, ax = plt.subplots(figsize=(5.5, 3))
        for c in (-2, -1, 0, 1, 2):
            ax.plot(xs, Fn(xs) + c, ":", color="0.6", lw=0.9)
        ax.plot(xs, Fn(xs) + C, color="navy", lw=2, label=f"C = {C}")
    ax.axhline(0, color="black", lw=0.6); ax.axvline(0, color="black", lw=0.6)
    ax.set_ylim(-10, 10); ax.set_xlabel("x"); ax.set_ylabel("F(x) + C"); ax.legend(); plt.show()


widgets.interact(calculadora_indefinida,
    f_txt=widgets.Text(value="2*x", description="f(x) ="),
    C=widgets.FloatSlider(value=0, min=-5, max=5, step=0.5, description="C"));

# %%
# Casos de prueba de la sección 1 (resultado conocido)
PRUEBAS_1 = [
    ("x^2 -> x^3/3", lambda: misma_funcion(antiderivada(x**2), x**3 / 3)),
    ("constante 3 -> 3x", lambda: misma_funcion(antiderivada(3 + 0 * x), 3 * x)),
    ("sqrt(x) -> (2/3) x^(3/2)", lambda: misma_funcion(antiderivada(sp.sqrt(x)), sp.Rational(2, 3) * x**sp.Rational(3, 2))),
    ("1/x -> ln|x| y aviso de la regla de la potencia",
     lambda: antiderivada(1 / x) == sp.log(sp.Abs(x)) and aviso_potencia(1 / x) != ""),
    ("entrada vacía o mal escrita se rechaza", lambda: _rechaza(lambda: parsear("")) and _rechaza(lambda: parsear("x^"))),
]
for nombre, prueba in PRUEBAS_1:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 1).** Calcula a mano $\int\big(4x^3-2x+1\big)dx$ y escribe tu respuesta **sin** la constante (la celda la agrega y compara derivando).

# %%
mi_F = None     # escribe tu antiderivada como texto, por ejemplo: "x**2 + 3*x"

if mi_F is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    f = 4 * x**3 - 2 * x + 1
    print("La calculadora da:", antiderivada(f), "+ C")
    try:
        ok = misma_funcion(sp.diff(parsear(str(mi_F)), x), f)
    except ValueError as err:
        print(err); ok = False
    print("coinciden (salvo la constante)" if ok else
          "NO coinciden: integra término a término con la regla de la potencia y deriva tu resultado para comprobar")
