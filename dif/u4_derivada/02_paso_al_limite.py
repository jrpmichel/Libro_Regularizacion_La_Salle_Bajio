# ID: DIF-U4-NB02
# Notebook: dif/u4_derivada.ipynb · sección 4.2 paso al límite
# Repositorio: dif/u4_derivada/02_paso_al_limite.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Paso al límite: definición de derivada
#
# Con $b=a+h$, la pendiente de la secante es el **cociente incremental** $\dfrac{f(a+h)-f(a)}{h}$. La derivada de $f$ en $a$ es su límite:
# $$f'(a)=\lim_{h\to 0}\frac{f(a+h)-f(a)}{h},$$
# si existe. Aquí $h$ es una **variable** que tiende a cero por los dos lados; un valor fijo de $h$ da solo la pendiente de una secante.
#
# La calculadora tabula el cociente con $h=\pm10^{-1},\dots,\pm10^{-k}$ y lo compara con el límite exacto que calcula `sympy`. Con $|h|$ menor que $10^{-8}$ el redondeo de la computadora empieza a dominar, por eso el control llega hasta $k=8$.

# %%
def tabla_cocientes(f, a, k=6):
    """[(h, cociente)] con h = ±10^-1, ..., ±10^-k, en punto flotante."""
    fn = sp.lambdify(x, f, "math")
    fa = evaluar(fn, float(a))
    if math.isnan(fa):
        raise ValueError(f"f no está definida en x = {a}: no hay punto de partida.")
    filas = []
    for j in range(1, k + 1):
        for hh in (10.0**-j, -10.0**-j):
            filas.append((hh, (evaluar(fn, float(a) + hh) - fa) / hh))
    return filas


def calculadora_limite(f_txt, a, k, n_cifras):
    try:
        f = parsear(f_txt)
        filas = tabla_cocientes(f, a, k)
        izq, der, d = derivada_en(f, a)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"{'h':>12}   cociente ({n_cifras} cifras significativas)")
    for hh, q in filas:
        print(f"{hh:>12g}   {cifras(q, n_cifras) if not math.isnan(q) else 'no definido'}")
    print("límite por la izquierda:", lado_texto(izq), "| por la derecha:", lado_texto(der))
    print("f'(a) =", a_texto(d) if d is not None else "no existe")
    print("Nota: k llega hasta 8 porque con |h| menor que 1e-8 el redondeo de la computadora domina el cociente.")


widgets.interact(calculadora_limite,
    f_txt=widgets.Text(value="20 + 70*exp(-0.1*x)", description="f(x) ="),
    a=widgets.FloatText(value=5.0, description="a"),
    k=widgets.IntSlider(value=4, min=1, max=8, description="k"),
    n_cifras=widgets.IntSlider(value=5, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 2 (resultado conocido)
cafe = 20 + 70 * sp.exp(-x / 10)
PRUEBAS_2 = [
    ("x² en 1: el cociente con h = 0.001 vale 2.001 y el límite es 2",
     lambda: cerca(dict(tabla_cocientes(x**2, 1, 3))[0.001], 2.001, 1e-9) and derivada_en(x**2, 1)[2] == 2),
    ("café en t = 5: límite -7e^(-1/2) ≈ -4.2457 °C/min",
     lambda: sp.simplify(derivada_en(cafe, 5)[2] + 7 * sp.exp(-sp.Rational(1, 2))) == 0),
    ("café: secante con h = 5 vale -3.3411", lambda: cerca(float(pendiente_secante(cafe, 5, 10)), -3.341117, 1e-5)),
    ("x³ en 2 con h = 0.1: 12.61", lambda: cerca(dict(tabla_cocientes(x**3, 2, 1))[0.1], 12.61, 1e-9)),
    ("1/x en 1 con h = 0.5: -2/3", lambda: pendiente_secante(1/x, 1, sp.Rational(3, 2)) == -sp.Rational(2, 3)),
    ("|x| en 0: laterales -1 y 1, la derivada no existe",
     lambda: derivada_en(sp.Abs(x), 0)[:2] == (-1, 1) and derivada_en(sp.Abs(x), 0)[2] is None),
    ("punto fuera del dominio se rechaza", lambda: _rechaza(lambda: tabla_cocientes(sp.sqrt(x), -1))),
    ("√x en 0: izquierda fuera del dominio, derecha +∞", lambda: derivada_en(sp.sqrt(x), 0) == (None, sp.oo, None)),
]
for nombre, prueba in PRUEBAS_2:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 2).** Calcula a mano el cociente incremental de $f(x)=x^2$ en $a=-2$ con $h=0.1$ y con $h=-0.1$. ¿A qué número se acercan? Compara con la calculadora.

# %%
mi_cociente_mas, mi_cociente_menos = None, None     # escribe dos números, por ejemplo: 1.5, 2.5

ref = dict(tabla_cocientes(x**2, -2, 1))
print("La calculadora da:", round(ref[0.1], 10), "y", round(ref[-0.1], 10), "| límite:", derivada_en(x**2, -2)[2])
print("falta tu cálculo" if mi_cociente_mas is None else
      ("coincide" if cerca(mi_cociente_mas, ref[0.1], 1e-6) and cerca(mi_cociente_menos, ref[-0.1], 1e-6)
       else "NO coincide: ((-2 + h)² - 4)/h = -4 + h"))
