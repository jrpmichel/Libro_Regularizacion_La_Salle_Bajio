# ID: DIF-U1-NB02
# Notebook: dif/u1_funciones.ipynb · sección 1.2 regla y gráfica
# Repositorio: dif/u1_funciones/02_regla_grafica.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Función como regla de asignación y como gráfica
#
# **Evaluar en una expresión.** Para $f(x)=x^2-3x+1$, el valor $f(a+1)$ se obtiene sustituyendo *cada* aparición de $x$ por $(a+1)$, con paréntesis. La primera calculadora lo hace con `sympy`.
#
# **Recta vertical.** La segunda calculadora cuenta en cuántos puntos corta la recta $x=c$ a una curva. Una sola recta con dos puntos basta para descartar que la curva sea gráfica de una función; una recta con un punto no basta para afirmar que lo es (hay que probar todas).

# %%
def evaluar_expresion(texto_f, texto_arg):
    f = parsear(texto_f)
    arg = parsear(texto_arg, variables=("x", "h", "a"))
    return sp.expand(f.subs(x, arg))


def calculadora_evaluar(texto_f, texto_arg):
    try:
        res = evaluar_expresion(texto_f, texto_arg)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"f(x) = {texto_f}")
    print(f"f({texto_arg}) = {res}")


widgets.interact(calculadora_evaluar,
    texto_f=widgets.Text(value="x**2 - 3*x + 1", description="f(x) ="),
    texto_arg=widgets.Text(value="a + 1", description="evaluar en"));

# %%
RELACIONES = {
    "x² + y² = 25 (circunferencia)": sp.Eq(x**2 + y**2, 25),
    "y² = x": sp.Eq(y**2, x),
    "y = x²": sp.Eq(y, x**2),
    "x² + 4y² = 16 (elipse)": sp.Eq(x**2 + 4*y**2, 16),
}


def recta_vertical(nombre, c):
    """Valores de y donde la recta x=c corta la relación (todos, con su signo)."""
    sols = sp.solve(RELACIONES[nombre].subs(x, sp.nsimplify(c)), y)
    return sorted({float(s) for s in sols if sp.im(s) == 0}, reverse=True)


def calculadora_recta_vertical(nombre, c, n_cifras):
    ys = recta_vertical(nombre, c)
    print(f"Relación: {nombre}      recta x = {cifras(c, n_cifras)}")
    print("Puntos de intersección:", [(cifras(c, n_cifras), cifras(v, n_cifras)) for v in ys] or "ninguno")
    if len(ys) > 1:
        print("Hay más de un punto: la curva NO es gráfica de y = f(x).")
    else:
        print("Con esta recta no hay violación. Eso no prueba que sea función: prueba otros valores de c.")
    t = np.linspace(-6, 6, 600)
    plt.figure(figsize=(4, 4)); plt.axvline(c, color="gray", ls="--")
    xs = np.linspace(-6, 6, 1200)
    for expr in sp.solve(RELACIONES[nombre], y):
        fn = sp.lambdify(x, expr, "numpy")
        with np.errstate(all="ignore"):
            yy = np.asarray(fn(xs), dtype=complex)
        yy = np.where(np.abs(yy.imag) < 1e-12, yy.real, np.nan)
        plt.plot(xs, yy, "b-")
    plt.plot([c]*len(ys), ys, "ro"); plt.xlim(-6, 6); plt.ylim(-6, 6); plt.grid(True)
    plt.gca().set_aspect("equal"); plt.xlabel("x"); plt.ylabel("y"); plt.show()


widgets.interact(calculadora_recta_vertical,
    nombre=widgets.Dropdown(options=list(RELACIONES), description="relación"),
    c=widgets.FloatSlider(value=3.0, min=-6, max=6, step=0.1, description="c"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
PRUEBAS_2 = [
    ("x²+y²=25 en x=3 da y=4 y y=-4",   lambda: recta_vertical("x² + y² = 25 (circunferencia)", 3) == [4.0, -4.0]),
    ("y²=x en x=0 da un solo punto",    lambda: recta_vertical("y² = x", 0) == [0.0]),
    ("y²=x en x=4 da y=2 y y=-2",       lambda: recta_vertical("y² = x", 4) == [2.0, -2.0]),
    ("y=x² en x=-1 da un solo punto",   lambda: recta_vertical("y = x²", -1) == [1.0]),
    ("f(a+1) para x²+2x-1 es a²+4a+2",  lambda: sp.simplify(evaluar_expresion("x**2 + 2*x - 1", "a + 1") - (a**2 + 4*a + 2)) == 0),
]
for nombre, prueba in PRUEBAS_2:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 2).** Calcula a mano $f(a+1)$ para $f(x)=x^2+2x-1$ y escríbelo en la celda. Anota: ¿coincide con la calculadora? Si no, ¿pusiste paréntesis al sustituir?

# %%
mi_resultado = ""      # por ejemplo: "a**2 + 4*a + 2"
comparar(mi_resultado, evaluar_expresion("x**2 + 2*x - 1", "a + 1"))
