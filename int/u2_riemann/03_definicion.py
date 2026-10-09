# ID: INT-U2-NB03
# Notebook: int/u2_riemann.ipynb · sección 2.3 la integral definida como límite de sumas
# Repositorio: int/u2_riemann/03_definicion.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. La integral definida como límite de sumas
#
# $$\int_a^b f(x)\,dx=\lim_{n\to\infty}\sum_{k=1}^{n}f(x_k^*)\,\Delta x,\qquad \Delta x=\frac{b-a}{n}.$$
#
# Para un polinomio, la suma derecha $R_n$ tiene fórmula cerrada en $n$ (con las fórmulas de $\sum k^p$) y el límite sale con álgebra. La calculadora escribe $R_n$, la simplifica, toma el límite y muestra $R_{10}$, $R_{100}$ y $R_{1000}$ para ver la convergencia en números.

# %%
def suma_derecha_simbolica(f, a, b):
    """R_n en forma cerrada (en n) y su límite, para f polinomial o simple."""
    a, b = sp.nsimplify(a), sp.nsimplify(b)
    if not a < b:
        raise ValueError("el intervalo debe cumplir a < b.")
    dx = (b - a) / n
    R = sp.simplify(sp.summation(f.subs(x, a + k * dx) * dx, (k, 1, n)))
    if R.has(sp.Sum):
        raise ValueError("sympy no encontró una fórmula cerrada para R_n; prueba con un polinomio.")
    return sp.factor(R), sp.limit(R, n, sp.oo)


def calculadora_limite(f_txt, a, b, n_cifras):
    try:
        f = parsear(f_txt)
        R, lim = suma_derecha_simbolica(f, a, b)
    except (ValueError, TypeError, NotImplementedError) as err:
        print("Revisa la entrada:", err); return
    an, bn = sp.nsimplify(a), sp.nsimplify(b)
    print("R_n =", sp.Sum(f.subs(x, an + k * (bn - an) / n) * (bn - an) / n, (k, 1, n)))
    print("    =", R)
    for m in (10, 100, 1000):
        print(f"  R_{m} ≈ {cifras(R.subs(n, m), n_cifras)}")
    print(f"límite cuando n → ∞: {lim} ≈ {cifras(lim, n_cifras)}   ({n_cifras} cifras significativas, redondeado)")


widgets.interact(calculadora_limite,
    f_txt=widgets.Text(value="x**2", description="f(x) ="),
    a=widgets.FloatText(value=0, description="a"), b=widgets.FloatText(value=1, description="b"),
    n_cifras=widgets.IntSlider(value=6, min=1, max=10, description="cifras sig."));

# %%
# Casos de prueba de la sección 3 (resultado conocido)
PRUEBAS_3 = [
    ("x^2 en [0, 1]: R_n = (n+1)(2n+1)/(6n^2) y límite 1/3",
     lambda: sp.simplify(suma_derecha_simbolica(x**2, 0, 1)[0] - (n + 1) * (2 * n + 1) / (6 * n**2)) == 0
             and suma_derecha_simbolica(x**2, 0, 1)[1] == sp.Rational(1, 3)),
    ("x en [0, 3]: límite 9/2 (triángulo)", lambda: suma_derecha_simbolica(x, 0, 3)[1] == sp.Rational(9, 2)),
    ("x^3 en [0, 2]: límite 4", lambda: suma_derecha_simbolica(x**3, 0, 2)[1] == 4),
    ("constante 5 en [1, 4]: 15", lambda: suma_derecha_simbolica(5 + 0 * x, 1, 4)[1] == 15),
    ("a >= b se rechaza", lambda: _rechaza(lambda: suma_derecha_simbolica(x, 2, 2))),
]
for nombre, prueba in PRUEBAS_3:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 3).** Escribe a mano $R_n$ para $f(x)=2x$ en $[0,1]$, simplifícala con la fórmula de $\sum k$ y toma el límite.

# %%
mi_Rn = None       # escribe la expresión en n como texto, por ejemplo: "3 - 1/n"
mi_limite = None   # escribe un número

if mi_Rn is None or mi_limite is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    R, lim = suma_derecha_simbolica(2 * x, 0, 1)
    print("La calculadora da: R_n =", R, "  límite =", lim)
    try:
        ok = sp.simplify(parsear(str(mi_Rn), variables=("n",)) - R) == 0 and cerca(mi_limite, lim, 1e-9)
    except ValueError as err:
        print(err); ok = False
    print("coinciden" if ok else "NO coinciden: R_n = suma de 2(k/n)(1/n) = (2/n^2)·n(n+1)/2")
