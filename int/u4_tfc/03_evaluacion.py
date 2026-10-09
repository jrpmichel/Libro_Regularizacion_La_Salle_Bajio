# ID: INT-U4-NB03
# Notebook: int/u4_tfc.ipynb · sección 4.3 evaluación de integrales definidas
# Repositorio: int/u4_tfc/03_evaluacion.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Integral con signo y área total
#
# $$\text{integral con signo}=\int_a^bf(x)\,dx,\qquad \text{área total}=\int_a^b|f(x)|\,dx=\sum_k\Big|\int_{c_{k}}^{c_{k+1}}f(x)\,dx\Big|,$$
#
# con $c_k$ los ceros de $f$ dentro de $(a,b)$. La calculadora encuentra los ceros, integra cada región con el teorema fundamental y entrega las dos cantidades. Decide tú cuál pide el problema.

# %%
def integral_y_area(f, a, b):
    """(integral con signo, área total, cortes, integrales por región)."""
    if not a < b:
        raise ValueError("el intervalo debe cumplir a < b.")
    cortes = [sp.nsimplify(a)] + ceros_en(f, a, b) + [sp.nsimplify(b)]
    partes = [evaluar(f, p, q)[1] for p, q in zip(cortes, cortes[1:])]
    return sp.nsimplify(sum(partes)), sp.nsimplify(sum(abs(v) for v in partes)), cortes, partes


def calculadora_area(f_txt, a, b, n_cifras):
    try:
        f = parsear(f_txt)
        integral, area, cortes, partes = integral_y_area(f, a, b)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print("ceros dentro del intervalo:", ", ".join(cifras(c, n_cifras) for c in cortes[1:-1]) or "ninguno")
    for p, q, v in zip(cortes, cortes[1:], partes):
        exacto = f"{v} ≈ " if len(str(v)) < 30 else ""
        print(f"  ∫ de {cifras(p, n_cifras)} a {cifras(q, n_cifras)} = {exacto}{cifras(v, n_cifras)}")
    print(f"integral con signo ≈ {cifras(integral, n_cifras)}   área total ≈ {cifras(area, n_cifras)}"
          f"   ({n_cifras} cifras significativas)")
    xs = np.linspace(a, b, 400); y = numerica(f)(xs)
    fig, ax = plt.subplots(figsize=(5.5, 3))
    ax.fill_between(xs, 0, y, where=y >= 0, color="0.8", label="cuenta positiva")
    ax.fill_between(xs, 0, y, where=y <= 0, facecolor="white", edgecolor="black", hatch="///", label="cuenta negativa")
    ax.plot(xs, y, color="navy"); ax.axhline(0, color="black", lw=0.6); ax.legend(); plt.show()


widgets.interact(calculadora_area,
    f_txt=widgets.Text(value="sin(x)", description="f(x) ="),
    a=widgets.FloatText(value=0, description="a"), b=widgets.FloatText(value=4.712389, description="b"),
    n_cifras=widgets.IntSlider(value=5, min=1, max=10, description="cifras sig."));

# %%
# Casos de prueba de la sección 3 (resultado conocido)
PRUEBAS_3 = [
    ("sen x en [0, 3π/2]: integral 1 y área 3",
     lambda: integral_y_area(sp.sin(x), 0, 3 * sp.pi / 2)[:2] == (1, 3)),
    ("x^3 en [-1, 2]: integral 15/4 y área 17/4",
     lambda: integral_y_area(x**3, -1, 2)[:2] == (sp.Rational(15, 4), sp.Rational(17, 4))),
    ("constante 2 en [0, 3]: 6 y 6", lambda: integral_y_area(2 + 0 * x, 0, 3)[:2] == (6, 6)),
    ("cos x en [0, 2π]: integral 0 y área 4", lambda: integral_y_area(sp.cos(x), 0, 2 * sp.pi)[:2] == (0, 4)),
    ("a >= b se rechaza", lambda: _rechaza(lambda: integral_y_area(x, 2, 1))),
]
for nombre, prueba in PRUEBAS_3:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 3).** Calcula a mano $\int_0^2(x-1)\,dx$ y el área total entre $y=x-1$ y el eje en $[0,2]$.

# %%
mi_integral = None     # escribe un número
mi_area = None         # escribe un número

if mi_integral is None or mi_area is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    integral, area, *_ = integral_y_area(x - 1, 0, 2)
    print("La calculadora da: integral =", integral, "  área =", area)
    print("coinciden" if cerca(mi_integral, integral) and cerca(mi_area, area) else
          "NO coinciden: los dos triángulos miden 1/2; con signo se restan y en el área se suman")
