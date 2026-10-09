# ID: DIF-U2-NB05
# Notebook: dif/u2_limite.ipynb · sección 2.5 sándwich y límites trigonométricos
# Repositorio: dif/u2_limite/05_sandwich_trig.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 5. Teorema del sándwich y límites trigonométricos
#
# El límite $\lim_{x\to0}\dfrac{\mathrm{sen}\,x}{x}=1$ vale **con $x$ en radianes**. Con $x$ en grados el mismo cociente tiende a $\pi/180\approx0.01745$. La primera calculadora tabula $\dfrac{\mathrm{sen}(kx)}{x}$ con la unidad del ángulo como control visible.
#
# La segunda dibuja una función atrapada entre dos cotas, que es lo que pide el teorema del sándwich: si $g(x)\le f(x)\le h(x)$ cerca de $a$ y $g$ y $h$ tienden al mismo $L$, entonces $f$ también.

# %%
def seno_kx_sobre_x(k, unidad, v):
    """sen(k·v)/v con v en la unidad indicada ('radianes' o 'grados')."""
    if unidad not in ("radianes", "grados"):
        raise ValueError("La unidad debe ser 'radianes' o 'grados'.")
    if v == 0:
        raise ValueError("En x = 0 el cociente no está definido; el límite se explora cerca de 0.")
    ang = k * v if unidad == "radianes" else math.radians(k * v)
    return math.sin(ang) / v


def calculadora_seno(k, unidad, n_cifras):
    try:
        seno_kx_sobre_x(k, unidad, 0.1)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    esperado = k if unidad == "radianes" else k * math.pi / 180
    print(f"sen({cifras(k, 4)}·x)/x con x en {unidad}   (valores redondeados a {n_cifras} cifras significativas)")
    for j in range(1, 6):
        v = 10.0**-j
        print(f"   x = ±{v:<8} {cifras(seno_kx_sobre_x(k, unidad, v), n_cifras):>14}  {cifras(seno_kx_sobre_x(k, unidad, -v), n_cifras):>14}")
    print(f"límite: {cifras(esperado, n_cifras)}" + ("" if unidad == "radianes" else "  (= k·π/180: la unidad cambia el resultado)"))


widgets.interact(calculadora_seno,
    k=widgets.FloatText(value=3.0, description="k ="),
    unidad=widgets.Dropdown(options=["radianes", "grados"], value="radianes", description="x en"),
    n_cifras=widgets.IntSlider(value=6, min=1, max=10, description="cifras sig."));

# %%
_SANDWICH = {
    "x·sen(1/x) entre -|x| y |x|": (lambda v: v*np.sin(1/v), lambda v: -np.abs(v), lambda v: np.abs(v)),
    "x²·cos(1/x) entre -x² y x²": (lambda v: v**2*np.cos(1/v), lambda v: -v**2, lambda v: v**2),
    "sen(x)/x entre cos x y 1": (lambda v: np.sin(v)/v, lambda v: np.cos(v), lambda v: np.ones_like(v)),
}


def sandwich(caso, ancho):
    if caso not in _SANDWICH:
        print("Revisa la entrada: elige un caso del menú."); return
    if not 0 < ancho <= 2:
        print("Revisa la entrada: el ancho debe estar entre 0 y 2."); return
    f, g, h = _SANDWICH[caso]
    v = np.linspace(-ancho, ancho, 4001); v = v[np.abs(v) > 1e-6]
    ok = np.all(g(v) <= f(v) + 1e-12) and np.all(f(v) <= h(v) + 1e-12)
    print(caso, "| ¿la función queda entre las cotas en todo el intervalo dibujado?", "sí" if ok else "NO")
    plt.figure(figsize=(5, 3)); plt.plot(v, f(v), "b-", lw=0.8, label="f")
    plt.plot(v, g(v), "k:", label="cota inferior"); plt.plot(v, h(v), "k--", label="cota superior")
    plt.grid(True); plt.legend(fontsize=8); plt.xlabel("x"); plt.show()


widgets.interact(sandwich, caso=widgets.Dropdown(options=list(_SANDWICH), description="caso"),
                 ancho=widgets.FloatSlider(value=0.5, min=0.05, max=2, step=0.05, description="ancho"));

# %%
# Casos de prueba de la sección 5 (resultado conocido)
PRUEBAS_5 = [
    ("sen x / x -> 1 (radianes)", lambda: sp.limit(sp.sin(x)/x, x, 0) == 1),
    ("sen(3x)/x -> 3; la tabla con x = 1e-5 da 3.0000", lambda: round(seno_kx_sobre_x(3, "radianes", 1e-5), 4) == 3.0),
    ("en grados, sen x / x -> π/180 = 0.017453", lambda: cerca(seno_kx_sobre_x(1, "grados", 1e-5), math.pi/180, 1e-9)),
    ("(1 - cos x)/x² -> 1/2 y (1 - cos x)/x -> 0",
     lambda: sp.limit((1 - sp.cos(x))/x**2, x, 0) == sp.Rational(1, 2) and sp.limit((1 - sp.cos(x))/x, x, 0) == 0),
    ("tan x / x -> 1", lambda: sp.limit(sp.tan(x)/x, x, 0) == 1),
    ("x² cos(1/x) -> 0 por sándwich", lambda: sp.limit(x**2*sp.cos(1/x), x, 0) == 0),
    ("sen(π/x) vale 0 en x = 1/n, pero no tiene límite en 0",
     lambda: all(abs(math.sin(math.pi*n)) < 1e-9 for n in (10, 100, 1000)) and isinstance(sp.limit(sp.sin(sp.pi/x), x, 0), sp.AccumBounds)),
    ("una unidad inválida se rechaza", lambda: _rechaza(lambda: seno_kx_sobre_x(1, "gradianes", 0.1))),
]
for nombre, prueba in PRUEBAS_5:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 5).** Calcula a mano $\displaystyle\lim_{x\to0}\frac{\mathrm{sen}\,3x}{x}$ escribiendo $\dfrac{\mathrm{sen}\,3x}{x}=3\cdot\dfrac{\mathrm{sen}\,3x}{3x}$. Escribe tu resultado y compáralo con la calculadora en radianes. Después cambia a grados y anota por qué el número cambia.

# %%
mi_limite = None        # por ejemplo: 3

ref = sp.limit(sp.sin(3*x)/x, x, 0)
print("La calculadora da:", ref)
print("falta tu cálculo" if mi_limite is None else ("coincide" if cerca(mi_limite, float(ref)) else "NO coincide: ¿multiplicaste y dividiste por 3?"))
