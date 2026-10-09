# ID: DIF-U5-NB04
# Notebook: dif/u5_reglas.ipynb · sección 5.4 derivadas de orden superior
# Repositorio: dif/u5_reglas/04_orden_superior.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 4. Derivadas de orden superior
#
# $f''=(f')'$, $f'''=(f'')'$ y, en general, $f^{(n)}$. Con $y=f(x)$: $\dfrac{d^2y}{dx^2},\ \dfrac{d^3y}{dx^3},\dots$ Las unidades son las de $y$ entre las de $x$ elevadas a la $n$.
#
# | Posición $s(t)$ | $s'$ | $s''$ | $s'''$ |
# |---|---|---|---|
# | m | velocidad, m/s | aceleración, m/s² | sacudida (*jerk*), m/s³ |
#
# En una viga (modelo de Euler-Bernoulli): deflexión $y$, giro $y'$, momento $M=EI\,y''$, cortante $V=EI\,y'''$ y carga $EI\,y''''=-w$.
#
# La calculadora tabula $f, f', \dots, f^{(n)}$ y sus valores en un punto, y avisa cuando la derivada se vuelve cero.

# %%
def derivadas_sucesivas(f, n):
    """[f, f', ..., f^(n)] como expresiones simplificadas."""
    if not (isinstance(n, int) and n >= 1):
        raise ValueError("el orden n debe ser un entero mayor o igual que 1.")
    lista = [f]
    for _ in range(n):
        lista.append(sp.simplify(sp.diff(lista[-1], x)))
    return lista


def calculadora_orden(f_txt, n, a, n_cifras):
    try:
        lista = derivadas_sucesivas(parsear(f_txt), int(n))
    except ValueError as err:
        print("Revisa la entrada:", err); return
    for k, d in enumerate(lista):
        nombre = "f" + ("'" * k if k <= 3 else f"^({k})")
        try:
            val = cifras(evaluar_en(d, x, a), n_cifras)
        except ValueError:
            val = "no definida"
        print(f"{nombre:>7}(x) = {str(d):<32} en x = {a}: {val}")
        if d == 0:
            print("        desde aquí todas las derivadas son cero."); break


widgets.interact(calculadora_orden,
    f_txt=widgets.Text(value="2*x**5 - x**3 + 4*x", description="f(x) ="),
    n=widgets.IntSlider(value=6, min=1, max=8, description="orden n"),
    a=widgets.FloatText(value=1.0, description="a"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 4 (resultado conocido)
PRUEBAS_4 = [
    ("2x^5 - x^3 + 4x: quinta 240 y sexta 0",
     lambda: derivadas_sucesivas(parsear("2*x**5 - x**3 + 4*x"), 6)[5:] == [240, 0]),
    ("1/x: tercera -6/x^4", lambda: iguales(derivadas_sucesivas(1/x, 3)[3], -6/x**4)),
    ("4.9t²: segunda 9.8 (aceleración)", lambda: derivadas_sucesivas(parsear("4.9*x**2"), 2)[2] == sp.Rational(49, 5)),
    ("ménsula w = 2, L = 3, EI = 2000: EI·y'''' = -2",
     lambda: 2000 * derivadas_sucesivas(parsear("-(2/48000)*(x**4 - 12*x**3 + 54*x**2)"), 4)[4] == -2),
    ("x^4 - 6x²: segunda en 1 vale 0", lambda: derivadas_sucesivas(x**4 - 6*x**2, 2)[2].subs(x, 1) == 0),
    ("orden 0 se rechaza", lambda: _rechaza(lambda: derivadas_sucesivas(x**2, 0))),
]
for nombre, prueba in PRUEBAS_4:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 4).** Calcula a mano $f''(1)$ para $f(x)=x^3-2x^2+5$ y escribe el número.

# %%
mi_valor = None          # escribe un número

ref = derivadas_sucesivas(x**3 - 2*x**2 + 5, 2)[2].subs(x, 1)
if mi_valor is None:                     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    print("La calculadora da:", ref)
    print("coincide" if mi_valor == ref else "NO coincide: deriva dos veces antes de sustituir x = 1")
