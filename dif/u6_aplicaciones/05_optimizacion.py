# ID: DIF-U6-NB05
# Notebook: dif/u6_aplicaciones.ipynb · sección 6.5 optimización
# Repositorio: dif/u6_aplicaciones/05_optimizacion.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 5. Optimización: criterios de primera y segunda derivada
#
# **Intervalo cerrado.** Si $f$ es continua en $[a,b]$, el máximo y el mínimo absolutos están entre los puntos críticos de $(a,b)$ y los extremos $a$, $b$: se evalúa $f$ en todos y se compara. **Segunda derivada:** $f'(c)=0$ y $f''(c)<0$ da máximo local; $f''(c)>0$, mínimo; $f''(c)=0$, no decide.
#
# **Retoma tu predicción de la semilla.** La celda siguiente optimiza la caja de cartón de 60 cm × 40 cm.

# %%
def extremos_intervalo(f, a, b):
    """(candidatos {x: f(x)}, x del máximo, x del mínimo) en [a, b]."""
    a, b = sp.nsimplify(a), sp.nsimplify(b)
    if not a < b:
        raise ValueError("el intervalo debe cumplir a < b.")
    crit = sp.solveset(sp.diff(f, x), x, domain=sp.Interval.open(a, b))
    if crit is sp.S.EmptySet:
        crit = []
    elif not isinstance(crit, sp.FiniteSet):
        raise ValueError("no pude listar los puntos críticos en el intervalo; prueba con otra función.")
    cand = {c: sp.simplify(f.subs(x, c)) for c in [a, b, *crit]}
    return cand, max(cand, key=lambda c: float(cand[c])), min(cand, key=lambda c: float(cand[c]))


def clasificar(f, c):
    """'máximo local', 'mínimo local' o 'sin extremo' en un punto crítico c (segunda derivada o cambio de signo)."""
    d2 = sp.diff(f, x, 2).subs(x, c)
    if d2 < 0:
        return "máximo local"
    if d2 > 0:
        return "mínimo local"
    d1 = sp.diff(f, x)
    e = sp.Rational(1, 1000)
    izq, der = d1.subs(x, c - e), d1.subs(x, c + e)
    return "máximo local" if izq > 0 > der else "mínimo local" if izq < 0 < der else "sin extremo"


caja = x*(60 - 2*x)*(40 - 2*x)
cand, xmax, _ = extremos_intervalo(caja, 0, 20)
print("caja: corte óptimo x ≈", cifras(xmax, 4), "cm, volumen ≈", cifras(cand[xmax], 4), "cm³")


def calculadora_optimizacion(f_txt, a, b, n_cifras):
    try:
        f = parsear(f_txt)
        cand, xM, xm = extremos_intervalo(f, a, b)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    for c, v in cand.items():
        nota = "extremo del intervalo" if c in (sp.nsimplify(a), sp.nsimplify(b)) else clasificar(f, c)
        print(f"  x = {cifras(c, n_cifras):>10}: f = {cifras(v, n_cifras):>10}   ({nota})")
    print(f"máximo absoluto en x ≈ {cifras(xM, n_cifras)}; mínimo absoluto en x ≈ {cifras(xm, n_cifras)}")


widgets.interact(calculadora_optimizacion,
    f_txt=widgets.Text(value="x**3 - 6*x**2 + 9*x + 1", description="f(x) ="),
    a=widgets.FloatText(value=0.5, description="a"), b=widgets.FloatText(value=4.5, description="b"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 5 (resultado conocido)
PRUEBAS_5 = [
    ("cúbica en [0.5, 4.5]: máximo en 4.5 (extremo), mínimo en 3",
     lambda: extremos_intervalo(x**3 - 6*x**2 + 9*x + 1, 0.5, 4.5)[1:] == (sp.Rational(9, 2), 3)),
    ("caja en [0, 20]: corte ≈ 7.85 cm, volumen ≈ 8450", lambda: cerca(xmax, 7.8475, 1e-4) and round(float(cand[xmax])) == 8450),
    ("caja con la máquina, [0, 5]: máximo en 5 con 7500", lambda: extremos_intervalo(caja, 0, 5)[1] == 5
                                                               and extremos_intervalo(caja, 0, 5)[0][5] == 7500),
    ("x(10 - x): cuadrado de 5, área 25", lambda: extremos_intervalo(x*(10 - x), 0, 10)[1] == 5),
    ("x^4 en 0: mínimo aunque f''(0) = 0", lambda: clasificar(x**4, 0) == "mínimo local"),
    ("x^3 en 0: sin extremo", lambda: clasificar(x**3, 0) == "sin extremo"),
    ("intervalo con a ≥ b se rechaza", lambda: _rechaza(lambda: extremos_intervalo(x**2, 3, 1))),
]
for nombre, prueba in PRUEBAS_5:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 5).** Calcula a mano el volumen de la caja con los cortes de la semilla, 2 cm, 8 cm y 15 cm, y compáralos con el óptimo.

# %%
mis_V = None             # escribe una lista de tres números: [V(2), V(8), V(15)]

if mis_V is None:                     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    print("La calculadora da:", [caja.subs(x, c) for c in (2, 8, 15)], "| óptimo ≈", cifras(cand[xmax], 4), "cm³")
    print("coinciden" if list(mis_V) == [caja.subs(x, c) for c in (2, 8, 15)] else "NO coinciden: V = x(60 - 2x)(40 - 2x)")
