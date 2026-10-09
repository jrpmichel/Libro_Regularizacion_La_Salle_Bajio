# ID: DIF-U5-NB01
# Notebook: dif/u5_reglas.ipynb · sección 5.1 constante, potencia y suma
# Repositorio: dif/u5_reglas/01_potencia.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Reglas de constante, potencia y suma
#
# | Regla | Fórmula | Condición |
# |---|---|---|
# | constante | $(c)'=0$ | |
# | potencia | $(x^n)'=n\,x^{n-1}$ | $n$ entero ($x\neq0$ si $n<0$) o $n$ real con $x>0$; exponente **constante** |
# | múltiplo constante | $(c\,f)'=c\,f'$ | |
# | suma y resta | $(f\pm g)'=f'\pm g'$ | |
#
# Antes de derivar, pasa raíces y cocientes a potencias: $\sqrt{x}=x^{1/2}$, $5/x^4=5x^{-4}$.
#
# La calculadora deriva una suma de potencias término por término y dice qué regla usó en cada uno. Si un término no es de la forma $c\,x^n$ con $n$ constante (por ejemplo $2^x$, $x^x$ o $x\,e^x$), lo rechaza: necesita otra regla.

# %%
def terminos_potencia(f):
    """Lista de (término, regla, derivada) de una suma de términos c*x**n con n constante."""
    filas = []
    for t in sp.Add.make_args(sp.expand(f)):
        if not t.has(x):
            filas.append((t, "constante", sp.Integer(0)))
            continue
        c, n = t.as_coeff_exponent(x)
        if c.has(x) or n.has(x):
            raise ValueError(f"el término {t} no es de la forma c·x^n con exponente constante: "
                             "necesita la regla del producto, de la cadena o la derivación logarítmica.")
        regla = "potencia" if c == 1 else "múltiplo constante y potencia"
        filas.append((t, regla, sp.simplify(c * n * x**(n - 1))))
    return filas


def calculadora_potencia(f_txt, a, n_cifras):
    try:
        f = parsear(f_txt)
        filas = terminos_potencia(f)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    for t, regla, d in filas:
        print(f"  {str(t):>16}  ->  {str(d):<18} ({regla})")
    total = sp.Add(*[d for _, _, d in filas])
    print("f'(x) =", total, "   (regla de la suma)")
    print("comprobación con sp.diff:", "coincide" if iguales(total, sp.diff(f, x)) else "NO coincide")
    try:
        va = evaluar_en(total, x, a)
        print(f"f'({a}) = {va} ≈ {cifras(va, n_cifras)}   | cociente con h = 1e-6: {cifras(cociente_num(f, x, a), n_cifras)}"
              f"   ({n_cifras} cifras significativas)")
    except (ValueError, ZeroDivisionError) as err:
        print("En a:", err)


widgets.interact(calculadora_potencia,
    f_txt=widgets.Text(value="4*x**5 - 3/x**2 + 7*sqrt(x) - 2", description="f(x) ="),
    a=widgets.FloatText(value=1.0, description="a"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 1 (resultado conocido)
def _derivada_terminos(texto):
    return sp.Add(*[d for _, _, d in terminos_potencia(parsear(texto))])

PRUEBAS_1 = [
    ("4x^5 - 3/x^2 + 7√x - 2 -> 20x^4 + 6/x^3 + 7/(2√x)",
     lambda: iguales(_derivada_terminos("4*x**5 - 3/x**2 + 7*sqrt(x) - 2"), 20*x**4 + 6/x**3 + sp.Rational(7, 2)/sp.sqrt(x))),
    ("x^3 - 3x -> 3x^2 - 3", lambda: iguales(_derivada_terminos("x**3 - 3*x"), 3*x**2 - 3)),
    ("constante 7 -> 0", lambda: _derivada_terminos("7") == 0),
    ("6√x - 2/x en 1: 5", lambda: _derivada_terminos("6*sqrt(x) - 2/x").subs(x, 1) == 5),
    ("2^x se rechaza (exponente variable)", lambda: _rechaza(lambda: terminos_potencia(parsear("2**x")))),
    ("x·e^x se rechaza (producto)", lambda: _rechaza(lambda: terminos_potencia(parsear("x*exp(x)")))),
]
for nombre, prueba in PRUEBAS_1:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 1).** Deriva a mano $f(x)=2x^4-\dfrac{5}{x}+3\sqrt{x}$. Escribe tu resultado como texto (con el formato de `"3*x**2 + 1/x"`) y compáralo con la calculadora.

# %%
mi_derivada = None       # escribe un texto entre comillas

contrasta(mi_derivada, sp.diff(parsear("2*x**4 - 5/x + 3*sqrt(x)"), x),
          pista="¿Pasaste 5/x a 5x^(-1) antes de derivar? ¿El signo del segundo término quedó positivo?")
