# ID: DIF-U2-NB01
# Notebook: dif/u2_limite.ipynb · sección 2.1 tablas desde ambos lados
# Repositorio: dif/u2_limite/01_tablas.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Aproximación numérica por tablas desde ambos lados
#
# Para explorar $\lim_{x\to a} f(x)$ se evalúa $f$ en $a-h$ y en $a+h$ con $h=0.1,\ 0.01,\ 0.001,\dots$ El valor $f(a)$ no interviene: muchas veces ni siquiera existe.
#
# La calculadora siguiente arma esa tabla. Si los dos lados se acercan al mismo número, la tabla **sugiere** ese límite. Con $h$ menor que $10^{-5}$ aparece otro riesgo, el redondeo de la máquina (Laboratorio del Error 2.2): la calculadora te avisa y no usa esos renglones para su lectura.

# %%
def tabla_lados(texto, a, k=4):
    """Filas (h, f(a-h), f(a+h)) con h = 10^-1 ... 10^-k. Valida las entradas."""
    f_expr = parsear(texto)
    a = float(numero(a))
    if not 1 <= int(k) <= 12:
        raise ValueError("Elige entre 1 y 12 renglones.")
    f = sp.lambdify(x, f_expr, "numpy")
    filas = []
    for j in range(1, int(k) + 1):
        h = 10.0**-j
        filas.append((h, evaluar(f, a - h), evaluar(f, a + h)))
    return filas


def lectura(filas, tol=1e-3):
    """Qué sugiere la tabla, usando solo los renglones con h >= 1e-5."""
    utiles = [r for r in filas if r[0] >= 1e-5 * 0.999] or filas[:1]
    _, izq, der = utiles[-1]
    if math.isnan(izq) or math.isnan(der):
        return "f no está definida de algún lado: revisa el dominio."
    if abs(izq) > 1e6 or abs(der) > 1e6:
        return "los valores crecen sin límite: posible límite infinito (subtema 2.2)."
    if abs(izq - der) < tol * max(1, abs(izq)):
        return f"los dos lados se acercan a ≈ {cifras((izq + der) / 2, 4)} (sugerencia, no demostración)."
    return f"los lados no coinciden (≈ {cifras(izq, 4)} y ≈ {cifras(der, 4)}): el límite no existe."


def calculadora_lados(texto, a, k, n_cifras):
    try:
        filas = tabla_lados(texto, a, k)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"f(x) = {texto}   cerca de a = {a}   (valores redondeados a {n_cifras} cifras significativas)")
    print(f"{'h':>8} | {'f(a - h)':>16} | {'f(a + h)':>16}")
    for h, izq, der in filas:
        print(f"{h:>8.0e} | {cifras(izq, n_cifras):>16} | {cifras(der, n_cifras):>16}")
    if filas[-1][0] < 1e-5 * 0.999:
        print("Aviso: con h < 1e-5 el redondeo de la máquina puede dominar (Laboratorio del Error 2.2);")
        print("       la lectura usa solo los renglones con h >= 1e-5.")
    print("Lectura:", lectura(filas))


widgets.interact(calculadora_lados,
    texto=widgets.Text(value="(x**3 - 1)/(x - 1)", description="f(x) ="),
    a=widgets.Text(value="1", description="a ="),
    k=widgets.IntSlider(value=4, min=1, max=12, description="renglones"),
    n_cifras=widgets.IntSlider(value=7, min=1, max=10, description="cifras sig."));

# %% [markdown]
# **La semilla, otra vez.** La rapidez media de la pelota entre $t=2$ y $t=2+h$ es $\dfrac{s(2+h)-s(2)}{h}$ con $s(t)=4.9t^2$. En la tabla, la variable se llama $x$ (es la $h$ de la semilla) y el punto es $a=0$. Compara con la predicción que anotaste.

# %%
calculadora_lados("(4.9*(2 + x)**2 - 4.9*4)/x", "0", 4, 6)

# %%
# Casos de prueba de la sección 1 (resultado conocido)
def ultimo_par(texto, a, k=4):
    return tabla_lados(texto, a, k)[-1][1:]

PRUEBAS_1 = [
    ("(x^3-1)/(x-1) cerca de 1: los dos lados tienden a 3",
     lambda: all(cerca(v, 3, 1e-3) for v in ultimo_par("(x**3 - 1)/(x - 1)", "1"))),
    ("f(0.9) = 2.71 y f(1.1) = 3.31",
     lambda: [round(v, 6) for v in tabla_lados("(x**3 - 1)/(x - 1)", "1", 1)[0][1:]] == [2.71, 3.31]),
    ("abs(x)/x cerca de 0: lados -1 y 1, la lectura dice que no existe",
     lambda: ultimo_par("abs(x)/x", "0") == (-1.0, 1.0) and "no existe" in lectura(tabla_lados("abs(x)/x", "0"))),
    ("semilla: la rapidez media tiende a 19.6 m/s",
     lambda: all(cerca(v, 19.6, 1e-3) for v in ultimo_par("(4.9*(2 + x)**2 - 4.9*4)/x", "0"))),
    ("factor P/A con n = 12 tiende a 12 cuando i -> 0",
     lambda: all(cerca(v, 12, 1e-2) for v in ultimo_par("(1 - (1 + x)**(-12))/x", "0"))),
    ("con h = 1e-9, (1 - cos x)/x^2 da 0 por redondeo (aviso del Laboratorio 2.2)",
     lambda: tabla_lados("(1 - cos(x))/x**2", "0", 9)[-1][1] == 0.0),
    ("una entrada con otra letra se rechaza con mensaje", lambda: _rechaza(lambda: tabla_lados("t**2", "0"))),
]
for nombre, prueba in PRUEBAS_1:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 1).** Calcula a mano $f(0.9)$ y $f(1.1)$ para $f(x)=\dfrac{x^3-1}{x-1}$ (sin calculadora: primero $0.9^3$). Escribe tus resultados y compáralos con la tabla. Anota: ¿cuál lado queda por debajo de 3 y cuál por arriba?

# %%
mi_f_09, mi_f_11 = None, None          # por ejemplo: 2.71, 3.31

esperado = tabla_lados("(x**3 - 1)/(x - 1)", "1", 1)[0][1:]
print("La calculadora da:", [round(v, 6) for v in esperado])
for nombre, mio, ref in zip(("f(0.9)", "f(1.1)"), (mi_f_09, mi_f_11), esperado):
    print(nombre, "->", "falta tu cálculo" if mio is None else ("coincide" if cerca(mio, ref, 1e-6) else "NO coincide"))
