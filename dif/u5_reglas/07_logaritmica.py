# ID: DIF-U5-NB07
# Notebook: dif/u5_reglas.ipynb · sección 5.7 derivación logarítmica
# Repositorio: dif/u5_reglas/07_logaritmica.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 7. Derivación logarítmica
#
# Si $y=f(x)>0$: $\ln y=\ln f(x)$, $\dfrac{y'}{y}=\big(\ln f(x)\big)'$ y $y'=f(x)\,\big(\ln f(x)\big)'$. Sirve para potencias con exponente variable, como $x^x$, y para productos y cocientes de muchas potencias.
#
# **Cambio relativo.** Si $y=k\,x_1^{a}x_2^{b}\cdots$, entonces $\dfrac{\Delta y}{y}\approx a\dfrac{\Delta x_1}{x_1}+b\dfrac{\Delta x_2}{x_2}+\cdots$ En el peor caso se suman los valores absolutos. El aporte de cada medición es su exponente por su error relativo.
#
# Dos calculadoras: la primera aplica la derivación logarítmica a una función positiva; la segunda reparte el error relativo de un resultado entre las mediciones.

# %%
def derivada_logaritmica(f_txt):
    """(ln f simplificado, (ln f)', f' = f·(ln f)')."""
    f = parsear(f_txt)
    xp = sp.symbols("xp", positive=True)        # el método supone f > 0; se trabaja con x > 0
    fp = f.subs(x, xp)
    lnf = sp.expand_log(sp.log(fp), force=True)
    dlnf = sp.simplify(sp.diff(lnf, xp))
    return lnf.subs(xp, x), dlnf.subs(xp, x), sp.simplify(fp * dlnf).subs(xp, x)


def sensibilidad(exponentes, errores_rel):
    """Aportes |a|·|e| de cada medición (ordenados de mayor a menor) y total en el peor caso."""
    if len(exponentes) != len(errores_rel) or not exponentes:
        raise ValueError("debe haber un exponente por cada error relativo.")
    aportes = sorted(((abs(a_) * abs(e_), i + 1) for i, (a_, e_) in enumerate(zip(exponentes, errores_rel))), reverse=True)
    return aportes, sum(a_ for a_, _ in aportes)


def calculadora_logaritmica(f_txt, a, n_cifras):
    try:
        lnf, dlnf, dfx = derivada_logaritmica(f_txt)
        va = evaluar_en(dfx, x, a)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print("ln f      =", lnf)
    print("(ln f)'   =", dlnf)
    print("f' = f·(ln f)' =", dfx)
    print(f"f'({a}) = {cifras(va, n_cifras)}  ({n_cifras} cifras significativas)"
          f"   | cociente con h = 1e-6: {cifras(cociente_num(parsear(f_txt), x, a), n_cifras)}")


widgets.interact(calculadora_logaritmica,
    f_txt=widgets.Text(value="x**x", description="f(x) ="),
    a=widgets.FloatText(value=2.0, description="a (> 0)"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));


def calculadora_sensibilidad(a1, e1, a2, e2, a3, e3):
    exps = [v for v in (a1, a2, a3)]
    errs = [v / 100 for v in (e1, e2, e3)]
    try:
        aportes, total = sensibilidad(exps, errs)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    for ap, i in aportes:
        print(f"medición {i}: exponente {exps[i - 1]:g} × error {100 * errs[i - 1]:.4g} % = {100 * ap:.4g} %")
    print(f"error relativo del resultado (peor caso): {100 * total:.4g} %")


print("Calculadora de sensibilidad: exponente y error relativo (%) de cada medición; usa 0 si no hay tercera.")
widgets.interact(calculadora_sensibilidad,
    a1=widgets.FloatText(value=2, description="exponente 1"), e1=widgets.FloatText(value=0.08, description="error 1 (%)"),
    a2=widgets.FloatText(value=1, description="exponente 2"), e2=widgets.FloatText(value=0.4167, description="error 2 (%)"),
    a3=widgets.FloatText(value=0, description="exponente 3"), e3=widgets.FloatText(value=0, description="error 3 (%)"));

# %%
# Casos de prueba de la sección 7 (resultado conocido)
PRUEBAS_7 = [
    ("x^x -> x^x(ln x + 1)", lambda: iguales(derivada_logaritmica("x**x")[2], x**x*(sp.log(x) + 1))),
    ("x^(ln x) en e: 2", lambda: cerca(evaluar_en(derivada_logaritmica("x**log(x)")[2], x, sp.E), 2, 1e-9)),
    ("(x+1)²√x/(x+3)³ en 1: 3/64",
     lambda: evaluar_en(derivada_logaritmica("(x + 1)**2*sqrt(x)/(x + 3)**3")[2], x, 1) == sp.Rational(3, 64)),
    ("eje actual: 2·0.08 % + 0.417 % ≈ 0.577 %", lambda: round(100*sensibilidad([2, 1], [0.0008, 0.5/120])[1], 3) == 0.577),
    ("con micrómetro: ≈ 0.433 %", lambda: round(100*sensibilidad([2, 1], [0.00008, 0.5/120])[1], 3) == 0.433),
    ("P = I²R con 2 % de corriente: 4 %", lambda: cerca(sensibilidad([2, 1], [0.02, 0])[1], 0.04)),
    ("listas de distinto tamaño se rechazan", lambda: _rechaza(lambda: sensibilidad([2, 1], [0.01]))),
]
for nombre, prueba in PRUEBAS_7:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 7).** Con $y=x^{\mathrm{sen}\,x}$, calcula a mano $y'(\pi/2)$ y escribe el número.

# %%
mi_valor = None          # escribe un número

ref = evaluar_en(derivada_logaritmica("x**sin(x)")[2], x, sp.pi/2)
if mi_valor is None:                     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    print("La calculadora da:", ref)
    print("coincide" if cerca(mi_valor, ref, 1e-6) else "NO coincide: ln y = sen x · ln x; deriva con la regla del producto")
