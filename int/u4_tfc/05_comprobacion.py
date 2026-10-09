# ID: INT-U4-NB05
# Notebook: int/u4_tfc.ipynb · sección 4.5 comprobación numérica contra la suma de Riemann
# Repositorio: int/u4_tfc/05_comprobacion.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 5. Comprobación numérica contra la suma de Riemann
#
# Si $F$ es correcta y $f$ continua, las diferencias entre $M_n$, $T_n$, $S_n$ y $F(b)-F(a)$ tienden a cero al crecer $n$. Si se **estancan**, $F$ está mal. Si las sumas **crecen sin cota**, la integral no existe como número finito.
#
# El comprobador recibe $f$ y, si quieres, tu antiderivada $F$ (si la dejas vacía usa la de `sympy`, o solo las sumas si no hay antiderivada elemental). Da la tabla para $n=2,4,\dots,256$ y un veredicto.

# %%
_NS = (2, 4, 8, 16, 32, 64, 128, 256)


def sumas(f, a, b, n):
    """(M_n, T_n, S_n) de f en [a, b]; n par."""
    g = numerica(f); xs = np.linspace(a, b, n + 1); y = g(xs); h = (b - a) / n
    M = float(g((xs[:-1] + xs[1:]) / 2).sum() * h)
    T = float(h * (y.sum() - (y[0] + y[-1]) / 2))
    S = float(h / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum()))
    return M, T, S


def comprobar(f, a, b, F=None):
    """(veredicto, valor de referencia o None, filas (n, M, T, S))."""
    if not a < b:
        raise ValueError("el intervalo debe cumplir a < b.")
    if not continua_en(f, a, b):
        g = numerica(f)
        Ms = []
        for n in (10, 100, 1000):
            xs = np.linspace(a, b, n + 1)
            with np.errstate(all="ignore"):
                Ms.append(float(g((xs[:-1] + xs[1:]) / 2).sum() * (b - a) / n))
        A = [abs(m) for m in Ms]
        crece = all(np.isfinite(Ms)) and A[2] > 3 * A[1] > 9 * A[0] > 0
        motivo = "f no es continua en [a, b] y las sumas crecen sin cota" if crece else "f no es continua en [a, b]"
        return f"no aplica: {motivo}; la integral no existe como número finito o el teorema no se puede usar", None, \
               [(n, m, float("nan"), float("nan")) for n, m in zip((10, 100, 1000), Ms)]
    filas = [(n,) + sumas(f, a, b, n) for n in _NS]
    if F is None:
        try:
            F = antiderivada(f)
        except ValueError:
            return "sin antiderivada elemental: el valor es el que dan las sumas", filas[-1][3], filas
    try:
        valor = valor_real(F, b) - valor_real(F, a)
    except (ValueError, TypeError) as err:
        raise ValueError("tu F no tiene valor real en los extremos del intervalo.") from err
    d = abs(filas[-1][3] - valor)                      # S_256 contra F(b) - F(a)
    e = abs(filas[-1][3] - filas[-2][3])               # cuánto cambió S_n en la última duplicación
    tol = 1e-9 * max(1, abs(valor))
    if d <= 10 * e + tol:
        return "coinciden", valor, filas
    if d > 100 * e + tol:
        return "no coinciden: revisa F", valor, filas
    return "no concluyente: las sumas aún no se estabilizan; aumenta n", valor, filas


def calculadora_comprobacion(f_txt, F_txt, a, b):
    try:
        f = parsear(f_txt)
        F = parsear(F_txt) if F_txt.strip() else None
        veredicto, valor, filas = comprobar(f, a, b, F)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    if valor is not None:
        print(f"valor de referencia: {cifras(valor, 8)}   (8 cifras significativas)")
    print("    n          M_n          T_n          S_n")
    for n, M, T, S in filas:
        celdas = [cifras(v, 8) if math.isfinite(v) else "—" for v in (M, T, S)]
        print(f"{n:5d}  {celdas[0]:>11}  {celdas[1]:>11}  {celdas[2]:>11}")
    print("Veredicto:", veredicto)


widgets.interact(calculadora_comprobacion,
    f_txt=widgets.Text(value="cos(2*x)", description="f(x) ="),
    F_txt=widgets.Text(value="sin(2*x)", description="tu F(x) ="),
    a=widgets.FloatText(value=0, description="a"), b=widgets.FloatText(value=1, description="b"));

# %%
# Casos de prueba de la sección 5 (resultado conocido)
PRUEBAS_5 = [
    ("cos 2x con F = sen(2x)/2: coinciden", lambda: comprobar(sp.cos(2 * x), 0, 1, sp.sin(2 * x) / 2)[0] == "coinciden"),
    ("cos 2x con F = sen 2x: no coinciden", lambda: comprobar(sp.cos(2 * x), 0, 1, sp.sin(2 * x))[0].startswith("no coinciden")),
    ("1/x^2 en [-1, 1]: no aplica, las sumas crecen", lambda: "crecen" in comprobar(1 / x**2, -1, 1)[0]),
    ("e^(-x^2) en [0, 1]: S_16 = 0.746824", lambda: round(sumas(sp.exp(-x**2), 0, 1, 16)[2], 6) == 0.746824),
    ("x^2 en [0, 3]: M_4 = 8.859375 (caso del prompt)", lambda: cerca(sumas(x**2, 0, 3, 4)[0], 8.859375)),
    ("sen 50x en [0, 10] con F correcta: coinciden aunque oscile",
     lambda: comprobar(sp.sin(50 * x), 0, 10, -sp.cos(50 * x) / 50)[0] == "coinciden"),
    ("√x en [0, 1] con F correcta: coinciden aunque S_n converja lento",
     lambda: comprobar(sp.sqrt(x), 0, 1, sp.Rational(2, 3) * x**sp.Rational(3, 2))[0].startswith("coinciden")),
]
for nombre, prueba in PRUEBAS_5:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 5).** Para $\int_0^1\big(x^3+1\big)dx$, calcula a mano $F(1)-F(0)$ y $M_2$.

# %%
mi_F = None    # escribe F(1) - F(0)
mi_M2 = None   # escribe M_2

if mi_F is None or mi_M2 is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    ref = evaluar(x**3 + 1, 0, 1)[1]; M2 = sumas(x**3 + 1, 0, 1, 2)[0]
    print("La calculadora da: F(1) - F(0) =", ref, "  M_2 =", M2)
    print("coinciden" if cerca(mi_F, ref) and cerca(mi_M2, M2) else
          "NO coinciden: M_2 evalúa x^3 + 1 en 0.25 y 0.75 y multiplica cada valor por 0.5")
