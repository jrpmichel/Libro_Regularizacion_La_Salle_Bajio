# ID: DIF-U3-NB03
# Notebook: dif/u3_transformaciones.ipynb · sección 3.3 estiramientos y compresiones
# Repositorio: dif/u3_transformaciones/03_estiramientos.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Estiramientos y compresiones
#
# - $a\,f(x)$ multiplica las alturas por $a$: estira si $|a|>1$, comprime si $0<|a|<1$.
# - $f(bx)$ divide las abscisas entre $b$: **comprime** si $|b|>1$ y estira si $0<|b|<1$.
#
# Para una función periódica de periodo $T$, $f(bx)$ tiene periodo $T/|b|$. La calculadora describe $y=a\,\mathrm{sen}\big(b(x-h)\big)+k$: amplitud, periodo, rango y primer máximo a partir de $x=h$.

# %%
def senoide(a, b, h=0, k=0):
    """Datos de y = a sen(b(x - h)) + k, en forma exacta con sympy."""
    a, b, h, k = (sp.nsimplify(v) for v in (a, b, h, k))
    if a == 0 or b == 0:
        raise ValueError("a y b deben ser distintos de cero (con a = 0 o b = 0 la función es constante).")
    T = 2 * sp.pi / abs(b)
    # primer máximo con x >= h: a·sen(u) es máximo en u = π/2 (a > 0) o en u = 3π/2 (a < 0),
    # con u = b(x - h); si b < 0 se toma el mismo ángulo un periodo antes para que x >= h
    u_max = sp.pi / 2 if a > 0 else 3 * sp.pi / 2
    xmax = h + (u_max if b > 0 else u_max - 2 * sp.pi) / b
    return {"amplitud": abs(a), "periodo": T, "rango": (k - abs(a), k + abs(a)), "primer máximo": sp.simplify(xmax)}


def calculadora_senoide(a, b, h, k, n_cifras):
    try:
        d = senoide(a, b, h, k)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"y = {cifras(a, 4)}·sen({cifras(b, 4)}(x - {cifras(h, 4)})) + {cifras(k, 4)}")
    print(f"amplitud = {cifras(d['amplitud'], n_cifras)}   periodo = {d['periodo']} ≈ {cifras(d['periodo'], n_cifras)}")
    print(f"rango = [{cifras(d['rango'][0], n_cifras)}, {cifras(d['rango'][1], n_cifras)}]   "
          f"primer máximo en x = {d['primer máximo']} ≈ {cifras(d['primer máximo'], n_cifras)}")
    xs = np.linspace(h, h + 2 * float(d["periodo"]), 600)
    plt.figure(figsize=(5, 3)); plt.plot(xs, a * np.sin(b * (xs - h)) + k, "b-")
    plt.axhline(k, color="gray", ls="--"); plt.grid(True); plt.xlabel("x (rad)"); plt.show()


widgets.interact(calculadora_senoide,
    a=widgets.FloatText(value=2.0, description="a"), b=widgets.FloatText(value=3.0, description="b"),
    h=widgets.FloatText(value=0.0, description="h"), k=widgets.FloatText(value=1.0, description="k"),
    n_cifras=widgets.IntSlider(value=5, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 3 (resultado conocido)
PRUEBAS_3 = [
    ("sen 2x tiene periodo π", lambda: senoide(1, 2)["periodo"] == sp.pi),
    ("2 sen 3x + 1: amplitud 2, periodo 2π/3, rango [-1, 3], primer máximo π/6",
     lambda: (lambda d: d["amplitud"] == 2 and d["periodo"] == 2*sp.pi/3 and d["rango"] == (-1, 3)
              and d["primer máximo"] == sp.pi/6)(senoide(2, 3, 0, 1))),
    ("sen(x/2) tiene periodo 4π", lambda: senoide(1, sp.Rational(1, 2))["periodo"] == 4*sp.pi),
    ("-3 sen x: amplitud 3, primer máximo en 3π/2", lambda: senoide(-3, 1)["primer máximo"] == 3*sp.pi/2),
    ("√(x/4) = √x/2 (estirar horizontalmente por 4 = comprimir verticalmente por 1/2)",
     lambda: sp.simplify(sp.sqrt(x/4) - sp.sqrt(x)/2) == 0),
    ("b = 0 se rechaza con mensaje", lambda: _rechaza(lambda: senoide(1, 0))),
]
for nombre, prueba in PRUEBAS_3:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 3).** Calcula a mano el periodo de $y=\mathrm{sen}(x/3)$ y di si la gráfica se comprimió o se estiró respecto a $\mathrm{sen}\,x$. Compara con la calculadora ($a=1$, $b=1/3$).

# %%
mi_periodo = None       # escribe un número, por ejemplo: 3.5 o 2*math.pi

ref = senoide(1, sp.Rational(1, 3))["periodo"]
print("La calculadora da:", ref, "≈", cifras(ref, 6))
print("falta tu cálculo" if mi_periodo is None else
      ("coincide" if cerca(mi_periodo, float(ref), 1e-6) else "NO coincide: el periodo de f(bx) es T/|b|"))
