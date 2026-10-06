# ID: DIF-U1-NB01
# Notebook: dif/u1_funciones.ipynb · sección 1.1 tablas
# Repositorio: dif/u1_funciones/01_tablas.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Qué es una función: tablas de valores $(x,y)$
#
# Una tabla define a $y$ como función de $x$ cuando **cada** valor de $x$ aparece con un solo valor de $y$. Dos $x$ distintos sí pueden compartir el mismo $y$.
#
# La calculadora siguiente construye la tabla de una regla $f$ en $n$ valores igualmente espaciados y calcula las diferencias primera y segunda. Si la primera diferencia es constante, los puntos están sobre una recta; si lo es la segunda, sobre una parábola (siempre que los datos sean exactos).

# %%
def tabla_valores(texto, x0, x1, n, n_cifras=4):
    """Devuelve (xs, ys) con ys=None donde f no está definida. Valida las entradas."""
    f = parsear(texto)
    if not (math.isfinite(x0) and math.isfinite(x1)) or x0 >= x1:
        raise ValueError("El valor inicial de x debe ser menor que el final.")
    if not 2 <= int(n) <= 50:
        raise ValueError("Elige entre 2 y 50 puntos.")
    fn = sp.lambdify(x, f, "numpy")
    xs = np.linspace(x0, x1, int(n))
    with np.errstate(all="ignore"):
        ys = fn(xs)
    ys = np.full_like(xs, float(ys)) if np.ndim(ys) == 0 else np.asarray(ys, dtype=float)
    return xs, np.where(np.isfinite(ys), ys, np.nan)


def calculadora_tabla(texto, x0, x1, n, n_cifras):
    try:
        xs, ys = tabla_valores(texto, x0, x1, n)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"f(x) = {texto}      (valores redondeados a {n_cifras} cifras significativas)")
    print(f"{'x':>10} | {'y':>12}")
    for u, v in zip(xs, ys):
        print(f"{cifras(u, n_cifras):>10} | {('no definida' if np.isnan(v) else cifras(v, n_cifras)):>12}")
    ok = ~np.isnan(ys)
    if ok.all() and len(ys) >= 4:
        d1, d2 = np.diff(ys), np.diff(ys, 2)
        print("\nprimera diferencia constante :", bool(np.allclose(d1, d1[0])))
        print("segunda diferencia constante :", bool(np.allclose(d2, d2[0])))
    plt.figure(figsize=(4.5, 3)); plt.plot(xs[ok], ys[ok], "o"); plt.grid(True)
    plt.xlabel("x"); plt.ylabel("y"); plt.title("Puntos de la tabla"); plt.show()


widgets.interact(calculadora_tabla,
    texto=widgets.Text(value="x**2 + 1", description="f(x) ="),
    x0=widgets.FloatText(value=0.0, description="x inicial"),
    x1=widgets.FloatText(value=4.0, description="x final"),
    n=widgets.IntSlider(value=5, min=2, max=30, description="puntos"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=10, description="cifras sig."));

# %%
def texto_a_pares(texto):
    """'0:20, 5:55, 5:60' -> [(0.0, 20.0), (5.0, 55.0), (5.0, 60.0)]"""
    pares = []
    for token in texto.replace(";", ",").split(","):
        if not token.strip():
            continue
        if token.count(":") != 1:
            raise ValueError(f"'{token.strip()}' no tiene el formato x:y (por ejemplo 5:55).")
        u, v = token.split(":")
        try:
            pares.append((float(u), float(v)))
        except ValueError:
            raise ValueError(f"'{token.strip()}' no es un par de números.") from None
    if not pares:
        raise ValueError("Escribe al menos un par x:y.")
    return pares


def es_funcion(pares):
    """(True, []) si cada x aparece con un solo y; si no, (False, pares en conflicto)."""
    vistos, conflictos = {}, []
    for u, v in pares:
        if u in vistos and vistos[u] != v:
            conflictos.append(((u, vistos[u]), (u, v)))
        vistos.setdefault(u, v)
    return len(conflictos) == 0, conflictos


def verificador_tabla(texto):
    try:
        pares = texto_a_pares(texto)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    for etiqueta, par in (("y es función de x", pares), ("x es función de y", [(v, u) for u, v in pares])):
        ok, conf = es_funcion(par)
        print(f"{etiqueta}: {'SÍ' if ok else 'NO'}", "" if ok else f"  (conflicto: {conf[0]})")
    print("Nota: con datos medidos, antes de concluir compara la diferencia entre lecturas "
          "con la incertidumbre del instrumento.")


widgets.interact(verificador_tabla, texto=widgets.Text(
    value="0:20, 5:55, 5:60, 10:88", description="pares x:y", layout=widgets.Layout(width="480px")));

# %%
# Casos de prueba de la sección 1 (resultado conocido)
def valor_en(texto, v):
    return float(parsear(texto).subs(x, v))

PRUEBAS_1 = [
    ("2*x + 1 en x=3 da 7",            lambda: valor_en("2*x + 1", 3) == 7),
    ("x**2 - 4 en x=-2 da 0",          lambda: valor_en("x**2 - 4", -2) == 0),
    ("3 - x en x=5 da -2",             lambda: valor_en("3 - x", 5) == -2),
    ("tabla repetida (1:4,2:6,2:7) no es función", lambda: not es_funcion(texto_a_pares("1:4, 2:6, 2:7"))[0]),
    ("la misma tabla sí define x como función de y", lambda: es_funcion([(v, u) for u, v in texto_a_pares("1:4, 2:6, 2:7")])[0]),
    ("y=x**2+1 en 0..4: segunda diferencia constante 2", lambda: np.allclose(np.diff(tabla_valores("x**2+1", 0, 4, 5)[1], 2), 2)),
]
for nombre, prueba in PRUEBAS_1:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 1).** Calcula a mano $f(-2)$, $f(0)$ y $f(3)$ para $f(x)=x^2-2x$ y compara con la tabla de la calculadora (usa $x$ de $-2$ a $3$ con 6 puntos). Anota: ¿coincidieron? Si no, ¿qué paso falló?

# %%
# Contrasta 1: escribe tus tres valores calculados a mano
mi_f_menos2, mi_f_0, mi_f_3 = None, None, None     # por ejemplo: 8, 0, 3

esperado = [valor_en("x**2 - 2*x", v) for v in (-2, 0, 3)]
print("La calculadora da:", esperado)
for nombre, mio, ok in zip(("f(-2)", "f(0)", "f(3)"), (mi_f_menos2, mi_f_0, mi_f_3), esperado):
    print(nombre, "->", "falta tu cálculo" if mio is None else ("coincide" if abs(mio - ok) < 1e-9 else "NO coincide"))
