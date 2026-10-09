# ID: DIF-U6-NB01
# Notebook: dif/u6_aplicaciones.ipynb · sección 6.1 crecimiento, concavidad e inflexión
# Repositorio: dif/u6_aplicaciones/01_grafica.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Análisis de la gráfica: crecimiento, concavidad y puntos de inflexión
#
# | Signo | Significado |
# |---|---|
# | $f'>0$ / $f'<0$ | $f$ crece / decrece |
# | $f''>0$ / $f''<0$ | cóncava hacia arriba / hacia abajo |
# | $f''$ cambia de signo en $c$ | punto de inflexión ($f'$ es máxima o mínima ahí) |
#
# $f'(c)=0$ no garantiza un extremo ($x^3$) y $f''(c)=0$ no garantiza inflexión ($x^4$): hay que revisar el **cambio de signo**.
#
# La calculadora encuentra los puntos críticos y los candidatos a inflexión (reales), revisa los signos a cada lado y arma la tabla de crecimiento y concavidad por intervalos.

# %%
def _signo(expr, v):
    val = float(expr.subs(x, v))
    return 1 if val > 1e-12 else (-1 if val < -1e-12 else 0)


def analizar(f, eps=sp.Rational(1, 1000)):
    """Diccionario con extremos locales [(c, tipo)] e inflexiones [c] de f (puntos donde f' o f'' se anulan)."""
    d1, d2 = sp.diff(f, x), sp.diff(f, x, 2)
    c1, c2 = sp.solveset(d1, x, domain=sp.S.Reals), sp.solveset(d2, x, domain=sp.S.Reals)
    if not all(isinstance(c, sp.FiniteSet) or c is sp.S.EmptySet for c in (c1, c2)):
        raise ValueError("no pude listar los puntos donde f' o f'' se anulan; prueba con un polinomio o una función más sencilla.")
    criticos, candidatos = sorted(c1, key=float), sorted(c2, key=float)
    extremos = []
    for c in criticos:
        izq, der = _signo(d1, c - eps), _signo(d1, c + eps)
        tipo = "máximo local" if (izq, der) == (1, -1) else "mínimo local" if (izq, der) == (-1, 1) else "sin extremo"
        extremos.append((c, tipo))
    inflexiones = [p for p in candidatos if _signo(d2, p - eps) * _signo(d2, p + eps) < 0]
    puntos = sorted(set(criticos) | set(candidatos), key=float)
    return {"f'": d1, "f''": d2, "extremos": extremos, "inflexiones": inflexiones, "puntos": puntos}


def tabla_signos(r):
    """[(a, b, signo de f', signo de f'')] en los intervalos que separan los puntos críticos y los candidatos a inflexión."""
    bordes = [-sp.oo, *r["puntos"], sp.oo]
    filas = []
    for a, b in zip(bordes[:-1], bordes[1:]):
        if a == -sp.oo and b == sp.oo:
            m = 0
        elif a == -sp.oo:
            m = b - 1
        elif b == sp.oo:
            m = a + 1
        else:
            m = (a + b) / 2
        filas.append((a, b, _signo(r["f'"], m), _signo(r["f''"], m)))
    return filas


def _txt(v, n):
    return "-∞" if v == -sp.oo else "∞" if v == sp.oo else cifras(v, n)


def calculadora_grafica(f_txt, n_cifras):
    try:
        f = parsear(f_txt)
        r = analizar(f)
    except (ValueError, TypeError) as err:
        print("Revisa la entrada:", err); return
    print("f'(x) =", sp.factor(r["f'"]), "   f''(x) =", sp.factor(r["f''"]))
    for c, tipo in r["extremos"]:
        print(f"  x = {c} ≈ {cifras(c, n_cifras)}: f'(x) = 0 -> {tipo}, f = {cifras(f.subs(x, c), n_cifras)}")
    for p in r["inflexiones"]:
        print(f"  x = {p} ≈ {cifras(p, n_cifras)}: f'' cambia de signo -> inflexión, f = {cifras(f.subs(x, p), n_cifras)}")
    if not r["inflexiones"]:
        print("  sin puntos de inflexión")
    print("tabla de signos:")
    for a, b, s1, s2 in tabla_signos(r):
        crece = "crece" if s1 > 0 else "decrece" if s1 < 0 else "constante"
        conc = "cóncava hacia arriba" if s2 > 0 else "cóncava hacia abajo" if s2 < 0 else "f'' = 0"
        print(f"  ({_txt(a, n_cifras)}, {_txt(b, n_cifras)}): f' {'+' if s1 > 0 else '-'} -> {crece};"
              f" f'' {'+' if s2 > 0 else '-'} -> {conc}")


widgets.interact(calculadora_grafica,
    f_txt=widgets.Text(value="x**3 - 6*x**2 + 9*x + 1", description="f(x) ="),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 1 (resultado conocido)
_cub = analizar(x**3 - 6*x**2 + 9*x + 1)
PRUEBAS_1 = [
    ("cúbica: máximo en 1 y mínimo en 3", lambda: _cub["extremos"] == [(1, "máximo local"), (3, "mínimo local")]),
    ("cúbica: inflexión en 2", lambda: _cub["inflexiones"] == [2]),
    ("x^4: mínimo en 0 y sin inflexión", lambda: analizar(x**4)["extremos"] == [(0, "mínimo local")]
                                                   and analizar(x**4)["inflexiones"] == []),
    ("x^3: punto crítico sin extremo", lambda: analizar(x**3)["extremos"] == [(0, "sin extremo")]),
    ("e^(-x²): inflexiones en ±1/√2", lambda: analizar(sp.exp(-x**2))["inflexiones"] == [-1/sp.sqrt(2), 1/sp.sqrt(2)]),
    ("horno 1 - 2e^-x + e^-2x: inflexión en ln 2", lambda: analizar(1 - 2*sp.exp(-x) + sp.exp(-2*x))["inflexiones"] == [sp.log(2)]),
    ("x^4 - 6x^2: inflexiones en -1 y 1; tabla de 6 intervalos", lambda: analizar(x**4 - 6*x**2)["inflexiones"] == [-1, 1]
                                                                    and len(tabla_signos(analizar(x**4 - 6*x**2))) == 6),
]
for nombre, prueba in PRUEBAS_1:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 1).** Encuentra a mano los puntos críticos y la inflexión de $f(x)=x^3-3x^2$. Escribe las listas.

# %%
mis_criticos = None      # escribe una lista de números, por ejemplo: [-1, 4]
mi_inflexion = None      # escribe un número

ref = analizar(x**3 - 3*x**2)
if mis_criticos is None:                     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    print("La calculadora da: críticos", [c for c, _ in ref["extremos"]], "| inflexión", ref["inflexiones"])
    print("coinciden" if sorted(mis_criticos) == [c for c, _ in ref["extremos"]] and [mi_inflexion] == ref["inflexiones"]
          else "NO coinciden: factoriza f' = 3x(x - 2) y resuelve f'' = 6x - 6 = 0")
