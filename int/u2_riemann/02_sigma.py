# ID: INT-U2-NB02
# Notebook: int/u2_riemann.ipynb · sección 2.2 notación sigma
# Repositorio: int/u2_riemann/02_sigma.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Notación sigma y fórmulas cerradas
#
# $$\sum_{k=1}^{n}k=\frac{n(n+1)}{2},\qquad \sum_{k=1}^{n}k^2=\frac{n(n+1)(2n+1)}{6},\qquad \sum_{k=1}^{n}k^3=\Big[\frac{n(n+1)}{2}\Big]^2,\qquad \sum_{k=1}^{n}c=n\,c.$$
#
# La calculadora recibe el término general en `k` y los límites. Si el límite superior es `n`, devuelve la fórmula cerrada y la **comprueba** con la suma término a término para $n=1,2,3$: una fórmula que falla con $n=1$ está mal, por elegante que se vea.

# %%
def suma_sigma(termino_txt, inferior, superior_txt):
    """Valor exacto de la suma (si el límite superior es un entero) o fórmula cerrada en n (si es 'n')."""
    t = parsear(termino_txt, variables=("k", "n"))
    if int(inferior) != inferior:
        raise ValueError("el límite inferior debe ser un entero.")
    inferior = int(inferior)
    sup = str(superior_txt).strip()
    if sup == "n":
        cerrada = sp.factor(sp.simplify(sp.summation(t, (k, inferior, n))))
        for m in range(max(1, inferior), max(1, inferior) + 3):        # comprobación término a término
            directa = sum(t.subs({k: j, n: m}) for j in range(inferior, m + 1))
            if sp.simplify(cerrada.subs(n, m) - directa) != 0:
                raise ValueError(f"la fórmula no coincide con la suma término a término en n = {m}.")
        return cerrada
    try:
        s_ = int(sup)
    except ValueError as err:
        raise ValueError("el límite superior debe ser un entero o la letra n.") from err
    if s_ < inferior:
        raise ValueError("el límite superior debe ser mayor o igual que el inferior.")
    if t.free_symbols - {k}:
        raise ValueError("con un límite superior numérico, el término solo puede depender de k.")
    return sp.nsimplify(sum(t.subs(k, j) for j in range(inferior, s_ + 1)))


def calculadora_sigma(termino_txt, inferior, superior_txt):
    try:
        r = suma_sigma(termino_txt, inferior, superior_txt)
    except (ValueError, TypeError) as err:
        print("Revisa la entrada:", err); return
    if str(superior_txt).strip() == "n":
        print(f"fórmula cerrada: {r}   (comprobada término a término para tres valores de n)")
    else:
        sup = int(str(superior_txt).strip())
        print(f"suma de {sup - int(inferior) + 1} términos = {r}")


widgets.interact(calculadora_sigma,
    termino_txt=widgets.Text(value="k**2", description="término"),
    inferior=widgets.IntText(value=1, description="desde k ="),
    superior_txt=widgets.Text(value="n", description="hasta"));

# %%
# Casos de prueba de la sección 2 (resultado conocido)
PRUEBAS_2 = [
    ("suma de k de 1 a 100 = 5050", lambda: suma_sigma("k", 1, "100") == 5050),
    ("suma de k^3 de 1 a n = (n(n+1)/2)^2", lambda: sp.simplify(suma_sigma("k**3", 1, "n") - (n * (n + 1) / 2)**2) == 0),
    ("suma de 3 de 1 a n = 3n (no 3)", lambda: suma_sigma("3", 1, "n") == 3 * n),
    ("suma de 3k^2 - 2k + 4 = n^3 + n^2/2 + 7n/2", lambda: sp.simplify(suma_sigma("3*k**2 - 2*k + 4", 1, "n")
                                                                       - (n**3 + n**2 / 2 + 7 * n / 2)) == 0),
    ("de k = 4 a 9 hay 6 términos; superior menor que inferior se rechaza",
     lambda: suma_sigma("1", 4, "9") == 6 and _rechaza(lambda: suma_sigma("k", 5, "2"))),
]
for nombre, prueba in PRUEBAS_2:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 2).** Calcula a mano $\displaystyle\sum_{k=1}^{6}\big(k^2-k\big)$ con las fórmulas cerradas.

# %%
mi_suma = None      # escribe un número

if mi_suma is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    ref = suma_sigma("k**2 - k", 1, "6")
    print("La calculadora da:", ref)
    print("coinciden" if mi_suma == ref else "NO coinciden: separa la suma en suma de k^2 menos suma de k, con n = 6")
