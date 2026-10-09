# ID: DIF-U2-NB02
# Notebook: dif/u2_limite.ipynb · sección 2.2 límites laterales y límites infinitos
# Repositorio: dif/u2_limite/02_laterales.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Lectura de límites en gráficas y límites laterales
#
# Una función por partes usa una regla a la izquierda de $a$ y otra a la derecha. El límite por la izquierda, $\lim_{x\to a^-}f(x)$, solo mira la regla de la izquierda; el límite por la derecha, $\lim_{x\to a^+}f(x)$, solo la de la derecha. El límite existe cuando los dos existen y son iguales. El valor $f(a)$ es un dato aparte.
#
# La calculadora recibe las dos reglas, el punto $a$ y el valor $f(a)$ (o la palabra `no`, si $f(a)$ no está definida), calcula los límites laterales con `sympy` y dibuja la función.

# %%
def laterales(izq_txt, der_txt, a):
    """(límite por la izquierda, límite por la derecha) de la función por partes en a."""
    fi, fd = parsear(izq_txt), parsear(der_txt)
    a = numero(a)
    return sp.limit(fi, x, a, "-"), sp.limit(fd, x, a, "+")


def veredicto(li, ld, valor=None):
    """Texto con la conclusión sobre el límite y la relación con f(a)."""
    if li.is_finite and ld.is_finite and sp.simplify(li - ld) == 0:
        txt = f"el límite existe y vale {a_texto(li)}."
        if valor is None:
            txt += " f(a) no está definida: hay un hueco."
        elif sp.simplify(valor - li) != 0:
            txt += f" Pero f(a) = {a_texto(valor)} es distinto del límite."
        else:
            txt += " Además coincide con f(a)."
        return txt
    if li == ld and li in (sp.oo, -sp.oo):
        return f"los dos lados tienden a {a_texto(li)}: el límite no existe como número (asíntota vertical)."
    return f"los laterales son distintos ({a_texto(li)} y {a_texto(ld)}): el límite no existe."


def calculadora_laterales(izq_txt, der_txt, a, valor_txt):
    try:
        li, ld = laterales(izq_txt, der_txt, a)
        valor = None if valor_txt.strip().lower() in ("no", "") else numero(valor_txt)
        a_num = float(numero(a))
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"x < {a}: f(x) = {izq_txt}      x > {a}: f(x) = {der_txt}      f({a}) = {valor_txt}")
    print("límite por la izquierda :", a_texto(li))
    print("límite por la derecha   :", a_texto(ld))
    print("conclusión              :", veredicto(li, ld, valor))
    fi = sp.lambdify(x, parsear(izq_txt), "numpy"); fd = sp.lambdify(x, parsear(der_txt), "numpy")
    xl = np.linspace(a_num - 3, a_num - 1e-3, 300); xr = np.linspace(a_num + 1e-3, a_num + 3, 300)
    with np.errstate(all="ignore"):
        yl = np.broadcast_to(fi(xl), xl.shape).astype(float); yr = np.broadcast_to(fd(xr), xr.shape).astype(float)
    yl = np.where(np.abs(yl) < 50, yl, np.nan); yr = np.where(np.abs(yr) < 50, yr, np.nan)
    plt.figure(figsize=(5, 3)); plt.plot(xl, yl, "b-", label="regla izquierda"); plt.plot(xr, yr, "b--", label="regla derecha")
    for li_, mk in ((li, "<"), (ld, ">")):
        if li_.is_finite:
            plt.plot([a_num], [float(li_)], "o", mfc="white", mec="b")
    if valor is not None:
        plt.plot([a_num], [float(valor)], "ro", label="f(a)")
    plt.axvline(a_num, color="gray", lw=0.6); plt.grid(True); plt.legend(fontsize=8)
    plt.xlabel("x"); plt.ylabel("f(x)"); plt.show()


widgets.interact(calculadora_laterales,
    izq_txt=widgets.Text(value="x + 3", description="izquierda"),
    der_txt=widgets.Text(value="2 + (x - 1)**2/9", description="derecha"),
    a=widgets.Text(value="-2", description="a ="),
    valor_txt=widgets.Text(value="3", description="f(a) ="));

# %%
# Casos de prueba de la sección 2 (resultado conocido)
def _lat(i, d, a):
    return laterales(i, d, a)

PRUEBAS_2 = [
    ("figura 2.2 en a = -2: izquierda 1, derecha 3", lambda: _lat("x + 3", "2 + (x - 1)**2/9", "-2") == (1, 3)),
    ("figura 2.2 en a = 1: el límite es 2 aunque g(1) = 4",
     lambda: "distinto" in veredicto(*_lat("2 + (x - 1)**2/9", "2 + (x - 1)**2/9", "1"), sp.Integer(4))),
    ("figura 2.2 en a = 4: izquierda 3, derecha +∞", lambda: _lat("2 + (x - 1)**2/9", "1/(x - 4)", "4") == (3, sp.oo)),
    ("abs(x)/x en 0: -1 y 1, no existe",
     lambda: _lat("abs(x)/x", "abs(x)/x", "0") == (-1, 1) and "no existe" in veredicto(-sp.Integer(1), sp.Integer(1))),
    ("cortante de la viga en x = 2: 8 y -4 kN, salto 12", lambda: (lambda l: l == (8, -4) and l[0] - l[1] == 12)(_lat("8", "-4", "2"))),
    ("1/(x - 2)**2 en 2: +∞ por ambos lados", lambda: _lat("1/(x - 2)**2", "1/(x - 2)**2", "2") == (sp.oo, sp.oo)),
    ("un número mal escrito se rechaza con mensaje", lambda: _rechaza(lambda: _lat("x", "x", "dos"))),
]
for nombre, prueba in PRUEBAS_2:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 2).** Para $f(x)=x+1$ si $x<2$ y $f(x)=x^2-1$ si $x\ge 2$, calcula a mano los dos límites laterales en $a=2$ y $f(2)$. Escríbelos y compáralos con la calculadora (`x + 1`, `x**2 - 1`, `2`, `3`). Anota si el límite existe y si coincide con $f(2)$.

# %%
mi_izq, mi_der, mi_valor = None, None, None      # por ejemplo: 3, 3, 3

li, ld = laterales("x + 1", "x**2 - 1", "2")
print("La calculadora da: izquierda", li, "| derecha", ld, "| f(2) =", 2**2 - 1)
for nombre, mio, ref in (("izquierda", mi_izq, li), ("derecha", mi_der, ld), ("f(2)", mi_valor, 3)):
    print(nombre, "->", "falta tu cálculo" if mio is None else ("coincide" if cerca(mio, float(ref)) else "NO coincide"))
