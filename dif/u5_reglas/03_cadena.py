# ID: DIF-U5-NB03
# Notebook: dif/u5_reglas.ipynb · sección 5.3 regla de la cadena
# Repositorio: dif/u5_reglas/03_cadena.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Regla de la cadena
#
# Si $y=f(u)$ y $u=g(x)$,
# $$\frac{dy}{dx}=\frac{dy}{du}\,\frac{du}{dx}=f'\big(g(x)\big)\,g'(x).$$
# Se deriva la función de afuera evaluada en la de adentro y se multiplica por la derivada de la de adentro.
#
# **Retoma tu predicción de la semilla.** La apuesta $10(x^2+1)^9$ da $5120$ en $x=1$. La celda siguiente calcula el cociente incremental con $h=0.001$ y la derivada con la regla de la cadena. Compara con lo que anotaste.

# %%
f_sem = (x**2 + 1)**10
print("cociente con h = 0.001:", round(cociente_num(f_sem, x, 1, 1e-3), 1))
print("apuesta 10(x²+1)^9    :", (10*(x**2 + 1)**9).subs(x, 1))
print("regla de la cadena    :", sp.diff(f_sem, x), "->", sp.diff(f_sem, x).subs(x, 1))

# %% [markdown]
# La calculadora recibe la función de afuera $f(u)$, escrita con la variable `u`, y la de adentro $g(x)$. Muestra cada paso, compara con `sp.diff` de la composición y evalúa en un punto junto al cociente incremental.

# %%
def cadena(f_txt, g_txt, a=None):
    """Pasos de la regla de la cadena para f(g(x)). Si se da a, también el valor en a."""
    F = parsear(f_txt, variables=("u",))
    G = parsear(g_txt, variables=("x",))
    dF, dG = sp.diff(F, u), sp.diff(G, x)
    resultado = sp.simplify(dF.subs(u, G) * dG)
    pasos = {"f'(u)": dF, "g'(x)": dG, "f'(g(x))": dF.subs(u, G), "(f(g(x)))'": resultado,
             "ok": iguales(resultado, sp.diff(F.subs(u, G), x))}
    if a is not None:
        ga = evaluar_en(G, x, a)
        try:
            evaluar_en(F, u, ga)
        except ValueError:
            raise ValueError(f"g({a}) = {ga} queda fuera del dominio de f.") from None
        pasos["valor"] = evaluar_en(resultado, x, a)
        pasos["cociente"] = cociente_num(F.subs(u, G), x, a)
    return pasos


def calculadora_cadena(f_txt, g_txt, a, n_cifras):
    try:
        p = cadena(f_txt, g_txt, a)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print("f'(u)       =", p["f'(u)"])
    print("g'(x)       =", p["g'(x)"])
    print("f'(g(x))    =", p["f'(g(x))"])
    print("derivada    =", p["(f(g(x)))'"], "  ->", "coincide con sp.diff" if p["ok"] else "NO coincide")
    print(f"en x = {a}: {cifras(p['valor'], n_cifras)}   | cociente con h = 1e-6: {cifras(p['cociente'], n_cifras)}"
          f"   ({n_cifras} cifras significativas; ángulos en radianes)")


widgets.interact(calculadora_cadena,
    f_txt=widgets.Text(value="sqrt(u)", description="f(u) ="),
    g_txt=widgets.Text(value="144*x**2 + 2500", description="g(x) ="),
    a=widgets.FloatText(value=10.0, description="a"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 3 (resultado conocido; los cuatro primeros son los del prompt del libro)
PRUEBAS_3 = [
    ("f = u^5, g = 3x + 1, a = 0: 15", lambda: cadena("u**5", "3*x + 1", 0)["valor"] == 15),
    ("f = √u, g = x² + 9, a = 4: 0.8", lambda: cadena("sqrt(u)", "x**2 + 9", 4)["valor"] == sp.Rational(4, 5)),
    ("f = u^10, g = x² + 1, a = 1: 10240", lambda: cadena("u**10", "x**2 + 1", 1)["valor"] == 10240),
    ("f = 1/u³, g = 2x - 1, a = 1: -6", lambda: cadena("1/u**3", "2*x - 1", 1)["valor"] == -6),
    ("dron: √u con 144x² + 2500 en 10: 144/13", lambda: cadena("sqrt(u)", "144*x**2 + 2500", 10)["valor"] == sp.Rational(144, 13)),
    ("g(a) fuera del dominio de f se rechaza", lambda: _rechaza(lambda: cadena("sqrt(u)", "-x**2 - 1", 0))),
    ("f escrita con x se rechaza", lambda: _rechaza(lambda: cadena("x**2", "3*x", 1))),
]
for nombre, prueba in PRUEBAS_3:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 3).** Deriva a mano $\sqrt{1+4x^2}$ y evalúa en $x=1$. Escribe el número con 4 cifras significativas.

# %%
mi_valor = None          # escribe un número

ref = cadena("sqrt(u)", "1 + 4*x**2", 1)["valor"]
if mi_valor is None:                     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    print("La calculadora da:", ref, "≈", cifras(ref, 4))
    print("coincide" if abs(mi_valor - float(ref)) < 5e-4 else "NO coincide: ¿multiplicaste por la derivada de 1 + 4x²?")
