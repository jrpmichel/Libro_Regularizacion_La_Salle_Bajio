# ID: DIF-U5-NB06
# Notebook: dif/u5_reglas.ipynb · sección 5.6 derivación implícita
# Repositorio: dif/u5_reglas/06_implicita.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6. Derivación implícita
#
# Para una curva dada por una ecuación en $x$ y $y$: se derivan los dos lados respecto a $x$ tratando $y$ como $y(x)$ (cada término con $y$ lleva un factor $y'$ por la cadena; $xy$ se deriva con la regla del producto) y se despeja $y'$. La pendiente depende del **punto**, no solo de $x$. Donde el coeficiente de $y'$ se anula, la fórmula no da la pendiente y puede haber una tangente vertical.
#
# La calculadora recibe la ecuación con un signo `=` y un punto. Comprueba que el punto esté en la curva y da $y'$, su valor y la recta tangente.

# %%
def implicita(ecuacion, px, py):
    """(y' como expresión, valor en el punto, (m, b) de la tangente y = m x + b)."""
    if ecuacion.count("=") != 1:
        raise ValueError("escribe una ecuación con un solo signo =, por ejemplo: x**2 + y**2 = 25")
    izq, der = ecuacion.split("=")
    F = parsear(izq, ("x", "y")) - parsear(der, ("x", "y"))
    punto = {x: sp.nsimplify(px), y: sp.nsimplify(py)}
    if sp.simplify(F.subs(punto)) != 0:
        raise ValueError(f"el punto ({px}, {py}) no está en la curva: el lado izquierdo menos el derecho vale "
                         f"{cifras(F.subs(punto), 4)}.")
    Fx, Fy = sp.diff(F, x), sp.diff(F, y)      # derivar F = 0 da Fx + Fy·y' = 0
    if sp.simplify(Fy.subs(punto)) == 0:
        raise ValueError("el coeficiente de y' se anula en ese punto: la fórmula no da la pendiente "
                         "(posible tangente vertical).")
    dydx = sp.simplify(-Fx / Fy)
    m = sp.simplify(dydx.subs(punto))
    return dydx, m, (m, sp.simplify(punto[y] - m * punto[x]))


def calculadora_implicita(ecuacion, px, py, n_cifras):
    try:
        dydx, m, (mm, b) = implicita(ecuacion, px, py)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print("y' =", dydx)
    print(f"en ({px}, {py}): y' = {m} ≈ {cifras(m, n_cifras)}")
    print(f"tangente: y = {cifras(mm, n_cifras)} x + {cifras(b, n_cifras)}   ({n_cifras} cifras significativas)")


widgets.interact(calculadora_implicita,
    ecuacion=widgets.Text(value="x**2 + y**2 = 25", description="ecuación"),
    px=widgets.FloatText(value=3.0, description="x0"), py=widgets.FloatText(value=4.0, description="y0"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 6 (resultado conocido)
PRUEBAS_6 = [
    ("circunferencia en (3, 4): -3/4", lambda: implicita("x**2 + y**2 = 25", 3, 4)[1] == sp.Rational(-3, 4)),
    ("circunferencia en (3, -4): 3/4", lambda: implicita("x**2 + y**2 = 25", 3, -4)[1] == sp.Rational(3, 4)),
    ("xy = 12 en (3, 4): -4/3", lambda: implicita("x*y = 12", 3, 4)[1] == sp.Rational(-4, 3)),
    ("folio x³ + y³ = 6xy en (3, 3): -1", lambda: implicita("x**3 + y**3 = 6*x*y", 3, 3)[1] == -1),
    ("x² + xy + y² = 7 en (1, 2): -4/5", lambda: implicita("x**2 + x*y + y**2 = 7", 1, 2)[1] == sp.Rational(-4, 5)),
    ("circunferencia en (5, 0): tangente vertical, se avisa", lambda: _rechaza(lambda: implicita("x**2 + y**2 = 25", 5, 0))),
    ("punto fuera de la curva se rechaza", lambda: _rechaza(lambda: implicita("x**2 + y**2 = 25", 1, 1))),
]
for nombre, prueba in PRUEBAS_6:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 6).** Calcula a mano la pendiente de $x^2+4y^2=8$ en $(2,1)$ y escríbela como número o fracción.

# %%
mi_pendiente = None      # escribe un número o un texto como "-2/3"

ref = implicita("x**2 + 4*y**2 = 8", 2, 1)[1]
if mi_pendiente is None:                 # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    print("La calculadora da:", ref)
    print("coincide" if iguales(sp.nsimplify(mi_pendiente), ref) else "NO coincide: d/dx(4y²) = 8y·y'")
