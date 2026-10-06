# ID: DIF-U1-NB05
# Notebook: dif/u1_funciones.ipynb · sección 1.5 lineal y cuadrática
# Repositorio: dif/u1_funciones/05_lineal_cuadratica.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 5. Funciones lineal y cuadrática
#
# Dos calculadoras. La primera obtiene pendiente y ordenada al origen de la recta que pasa por dos puntos. La segunda analiza $f(x)=ax^2+bx+c$: discriminante, **ambas** raíces con su signo, vértice, concavidad y rango. Si $a=0$ avisa que la función ya no es cuadrática.

# %%
def recta_dos_puntos(x1, y1, x2, y2):
    if x1 == x2:
        raise ValueError("Con x1 = x2 la recta es vertical y no es gráfica de una función de x.")
    m = (y2 - y1) / (x2 - x1)
    return m, y1 - m * x1


def calculadora_recta(x1, y1, x2, y2, n_cifras):
    try:
        m, b = recta_dos_puntos(x1, y1, x2, y2)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"pendiente m = {cifras(m, n_cifras)}      ordenada al origen b = {cifras(b, n_cifras)}")
    print(f"y = {cifras(m, n_cifras)} x + {cifras(b, n_cifras)}")


widgets.interact(calculadora_recta,
    x1=widgets.FloatText(value=1.0, description="x1"), y1=widgets.FloatText(value=3.0, description="y1"),
    x2=widgets.FloatText(value=4.0, description="x2"), y2=widgets.FloatText(value=12.0, description="y2"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
def cuadratica(A, B, C):
    """Análisis de a x² + b x + c. Devuelve un diccionario con todo lo que pide la teoría."""
    if A == 0:
        raise ValueError("a = 0: la función ya no es cuadrática, es la recta y = bx + c.")
    disc = B**2 - 4*A*C
    xv = -B / (2*A); yv = A*xv**2 + B*xv + C
    if disc > 0:
        raices = sorted([(-B + s*math.sqrt(disc)) / (2*A) for s in (1, -1)])
    elif disc == 0:
        raices = [xv]
    else:
        raices = []
    rango = f"[{cifras(yv, 6)}, ∞)" if A > 0 else f"(-∞, {cifras(yv, 6)}]"
    return {"disc": disc, "vertice": (xv, yv), "raices": raices, "rango": rango,
            "abre": "hacia arriba" if A > 0 else "hacia abajo", "f(0)": C}


def calculadora_cuadratica(A, B, C, n_cifras):
    try:
        r = cuadratica(A, B, C)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    xv, yv = r["vertice"]
    poli = sp.nsimplify(A)*x**2 + sp.nsimplify(B)*x + sp.nsimplify(C)
    print(f"f(x) = {poli}      abre {r['abre']}      f(0) = {cifras(r['f(0)'], n_cifras)}")
    print("discriminante:", cifras(r["disc"], n_cifras))
    if r["disc"] > 0:
        print("raíces (ambas):", ", ".join(cifras(v, n_cifras) for v in r["raices"]))
    elif r["disc"] == 0:
        print("una raíz real (doble):", cifras(r["raices"][0], n_cifras))
    else:
        print("sin raíces reales (el discriminante es negativo)")
    print("vértice:", (cifras(xv, n_cifras), cifras(yv, n_cifras)), "      rango:", r["rango"])
    ancho = max(2.0, abs(xv) + 3, *(abs(v - xv) + 1 for v in r["raices"]))
    xs = np.linspace(xv - ancho, xv + ancho, 400)
    plt.figure(figsize=(4.5, 3)); plt.plot(xs, A*xs**2 + B*xs + C, "b-")
    plt.plot([xv], [yv], "rs"); plt.plot(r["raices"], [0]*len(r["raices"]), "ko", mfc="white")
    plt.axvline(xv, color="gray", ls="--"); plt.axhline(0, color="k", lw=0.8); plt.grid(True)
    plt.xlabel("x"); plt.ylabel("y"); plt.show()


widgets.interact(calculadora_cuadratica,
    A=widgets.FloatText(value=1.0, description="a"), B=widgets.FloatText(value=-4.0, description="b"),
    C=widgets.FloatText(value=3.0, description="c"),
    n_cifras=widgets.IntSlider(value=5, min=1, max=8, description="cifras sig."));

# %%
def cerca(u, v, tol=1e-9):
    return abs(u - v) < tol

PRUEBAS_5 = [
    ("x²-4x+3: raíces 1 y 3, vértice (2,-1)",
     lambda: [cerca(p, q) for p, q in zip(cuadratica(1, -4, 3)["raices"], (1, 3))] == [True, True]
             and cuadratica(1, -4, 3)["vertice"] == (2.0, -1.0)),
    ("-2x²+400x-8000: raíces 22.54 y 177.46, vértice (100, 12000)",
     lambda: [round(v, 2) for v in cuadratica(-2, 400, -8000)["raices"]] == [22.54, 177.46]
             and cuadratica(-2, 400, -8000)["vertice"] == (100.0, 12000.0)),
    ("x²+2x+5: discriminante -16, sin raíces reales",
     lambda: cuadratica(1, 2, 5)["disc"] == -16 and cuadratica(1, 2, 5)["raices"] == []),
    ("x²-6x+9: una raíz doble en 3",   lambda: cuadratica(1, -6, 9)["raices"] == [3.0]),
    ("recta por (1,3) y (4,12): m=3, b=0", lambda: recta_dos_puntos(1, 3, 4, 12) == (3.0, 0.0)),
    ("a=0 se rechaza con mensaje",     lambda: _rechaza(lambda: cuadratica(0, 1, 1))),
    ("x1=x2 se rechaza con mensaje",   lambda: _rechaza(lambda: recta_dos_puntos(2, 1, 2, 5))),
]

def _rechaza(fn):
    try:
        fn()
    except ValueError:
        return True
    return False

for nombre, prueba in PRUEBAS_5:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 5).** Calcula a mano el discriminante, las dos raíces y el vértice de $f(x)=x^2-6x+5$. Compara con la calculadora. Anota: ¿escribiste las dos raíces? ¿el vértice queda en el punto medio de ellas?

# %%
mis_raices = []        # por ejemplo: [1, 5]
mi_vertice = None      # por ejemplo: (3, -4)

r = cuadratica(1, -6, 5)
print("La calculadora da: raíces", r["raices"], "| vértice", r["vertice"])
print("Tus raíces  :", mis_raices or "(falta tu cálculo a mano)")
print("Tu vértice  :", mi_vertice or "(falta tu cálculo a mano)")
if mis_raices:
    print("Raíces:", "coinciden" if sorted(mis_raices) == r["raices"] else "NO coinciden")

# %% [markdown]
# **Retoma tu predicción de la semilla.** Con los datos $(0\,^\circ\mathrm{C},\,0.50\ \mathrm{V})$ y $(30\,^\circ\mathrm{C},\,0.80\ \mathrm{V})$ de la tabla del sensor, obtén la recta $V(T)$ y calcula $V(25)$ y $V(100)$. Compara con lo que anotaste. La primera es una interpolación (dentro del intervalo medido); la segunda, una extrapolación, y solo vale si el sensor sigue siendo lineal fuera del rango calibrado.

# %%
m, b = recta_dos_puntos(0, 0.50, 30, 0.80)
print(f"V(T) = {m:.4g} T + {b:.4g}   [V, °C]")
print("V(25)  =", round(m*25 + b, 4), "V   (interpolación, dentro de 0 a 30 °C)")
print("V(100) =", round(m*100 + b, 4), "V  (extrapolación: supone que la recta sigue valiendo)")
