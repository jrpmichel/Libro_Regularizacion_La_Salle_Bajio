# ID: DIF-U5-NB05
# Notebook: dif/u5_reglas.ipynb · sección 5.5 exponencial, logaritmo y trigonométricas
# Repositorio: dif/u5_reglas/05_exp_log_trig.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 5. Derivadas de exponencial, logaritmo y funciones trigonométricas
#
# Con $x$ en **radianes** en las trigonométricas:
#
# | $f$ | $e^x$ | $a^x$ | $\ln x$ | $\mathrm{sen}\,x$ | $\cos x$ | $\tan x$ |
# |---|---|---|---|---|---|---|
# | $f'$ | $e^x$ | $a^x\ln a$ | $1/x$ ($x>0$) | $\cos x$ | $-\mathrm{sen}\,x$ | $1/\cos^2x$ |
#
# Con grados, $\mathrm{sen}(\theta^\circ)=\mathrm{sen}(\pi\theta/180)$ y la cadena agrega el factor $\pi/180\approx0.01745$.
#
# La calculadora tiene un selector visible de unidad angular. En grados reescribe cada `sin`, `cos` y `tan` con el argumento multiplicado por $\pi/180$ y muestra el resultado. Rechaza `log` de valores no positivos.

# %%
_GRADOS = {**_FUNCIONES,
           "sin": lambda e: sp.sin(sp.pi * e / 180), "cos": lambda e: sp.cos(sp.pi * e / 180),
           "tan": lambda e: sp.tan(sp.pi * e / 180)}


def derivar_angulos(texto, modo="radianes"):
    """(f, f') con los ángulos en el modo indicado ('radianes' o 'grados')."""
    if modo not in ("radianes", "grados"):
        raise ValueError("el modo debe ser 'radianes' o 'grados'.")
    f = parsear(texto, funciones=_GRADOS if modo == "grados" else None)
    return f, sp.diff(f, x)


def derivada_en(texto, a, modo="radianes"):
    """f'(a); ValueError si a queda fuera del dominio (por ejemplo, log de un valor no positivo)."""
    f, d = derivar_angulos(texto, modo)
    for arg in [g.args[0] for g in f.atoms(sp.log)]:
        if not evaluar_en(arg, x, a) > 0:
            raise ValueError(f"en x = {a} el argumento del logaritmo, {arg}, no es positivo.")
    evaluar_en(f, x, a)
    return evaluar_en(d, x, a)


def calculadora_exp_trig(f_txt, modo, a, n_cifras):
    try:
        f, d = derivar_angulos(f_txt, modo)
        va = derivada_en(f_txt, a, modo)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    unidad = ("por grado" if modo == "grados" else "por radián") if f.has(sp.sin, sp.cos, sp.tan) else ""
    print(f"modo angular: {modo.upper()}")
    print("f(x)  =", f)
    print("f'(x) =", sp.simplify(d))
    print(f"f'({a}) = {cifras(va, n_cifras)} {unidad}  ({n_cifras} cifras significativas)"
          f"   | cociente con h = 1e-6: {cifras(cociente_num(f, x, a), n_cifras)}")


widgets.interact(calculadora_exp_trig,
    f_txt=widgets.Text(value="sin(x)", description="f(x) ="),
    modo=widgets.ToggleButtons(options=["radianes", "grados"], description="ángulos"),
    a=widgets.FloatText(value=60.0, description="a"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 5 (resultado conocido)
PRUEBAS_5 = [
    ("e^x -> e^x", lambda: derivar_angulos("exp(x)")[1] == sp.exp(x)),
    ("2^x en 0: ln 2", lambda: iguales(derivada_en("2**x", 0), sp.log(2))),
    ("sen x en π/3 (radianes): 1/2", lambda: derivada_en("sin(x)", sp.pi/3) == sp.Rational(1, 2)),
    ("sen x en 60 (grados): π/360 ≈ 0.008727", lambda: iguales(derivada_en("sin(x)", 60, "grados"), sp.pi/360)),
    ("café 20 + 70e^(-0.1x) en 5: -4.2457", lambda: cerca(derivada_en("20 + 70*exp(-0.1*x)", 5), -4.2457146, 1e-6)),
    ("sigmoide en 0: 0.25", lambda: derivada_en("1/(1 + exp(-x))", 0) == sp.Rational(1, 4)),
    ("ln x en -1 se rechaza", lambda: _rechaza(lambda: derivada_en("log(x)", -1))),
    ("e^(-x) cos 2x en 0: -1", lambda: derivada_en("exp(-x)*cos(2*x)", 0) == -1),
]
for nombre, prueba in PRUEBAS_5:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 5).** Deriva a mano $e^{-2x}\,\mathrm{sen}\,3x$, evalúa en $x=0$ y escribe el número. Después cambia la calculadora a grados con la misma función y $a=0$: ¿por qué cambia el resultado?

# %%
mi_valor = None          # escribe un número

ref = derivada_en("exp(-2*x)*sin(3*x)", 0)
if mi_valor is None:                     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    print("La calculadora (radianes) da:", ref, "| en grados:", cifras(derivada_en("exp(-2*x)*sin(3*x)", 0, "grados"), 4))
    print("coincide" if mi_valor == ref else "NO coincide: regla del producto y cadena con la interior 3x")
