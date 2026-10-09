# ID: INT-U2-NB05
# Notebook: int/u2_riemann.ipynb · sección 2.5 trapecio, Simpson y tabla de error
# Repositorio: int/u2_riemann/05_trapecio_simpson.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 5. Regla del trapecio y regla de Simpson, con tabla de error contra $n$
#
# $$T_n=\frac{\Delta x}{2}\big[y_0+2y_1+\cdots+2y_{n-1}+y_n\big],\qquad S_n=\frac{\Delta x}{3}\big[y_0+4y_1+2y_2+\cdots+4y_{n-1}+y_n\big]\ \ (n\ \text{par}).$$
#
# Cotas: $|E_T|\le\dfrac{K_2(b-a)^3}{12n^2}$ y $|E_S|\le\dfrac{K_4(b-a)^5}{180n^4}$. Al duplicar $n$, el error del trapecio se divide entre $\approx4$ y el de Simpson entre $\approx16$.
#
# La calculadora acepta una **función** o una **tabla** de lecturas igualmente espaciadas. Con una tabla de $N$ lecturas hay $N-1$ franjas; Simpson exige que ese número sea par.

# %%
def trapecio_datos(y, h):
    y = np.asarray(y, dtype=float)
    if not h > 0:
        raise ValueError("el ancho h de cada franja debe ser positivo.")
    if y.size < 2:
        raise ValueError("se necesitan al menos dos lecturas.")
    return float(h * (y.sum() - (y[0] + y[-1]) / 2))


def simpson_datos(y, h):
    y = np.asarray(y, dtype=float)
    if not h > 0:
        raise ValueError("el ancho h de cada franja debe ser positivo.")
    franjas = y.size - 1
    if franjas < 2 or franjas % 2:
        raise ValueError(f"hay {franjas} franjas: Simpson necesita un número par (un número impar de lecturas). "
                         "Usa el trapecio o combina Simpson con un trapecio en la última franja.")
    return float(h / 3 * (y[0] + y[-1] + 4 * y[1:-1:2].sum() + 2 * y[2:-1:2].sum()))


def reglas_funcion(f, a, b, m):
    """(T_m, S_m) de la expresión f; S_m es None si m es impar."""
    if not a < b:
        raise ValueError("el intervalo debe cumplir a < b.")
    if int(m) != m or m < 1:
        raise ValueError("el número de franjas debe ser un entero positivo.")
    m = int(m)
    xs = np.linspace(a, b, m + 1)
    y = [valor_real(f, c) for c in xs]
    return trapecio_datos(y, (b - a) / m), (simpson_datos(y, (b - a) / m) if m % 2 == 0 else None)


def cota_maxima(expr, a, b):
    """Máximo de |expr| en [a, b], muestreado en 2001 puntos (para K_2 y K_4); inf si no está acotada."""
    with np.errstate(all="ignore"):
        v = np.abs(numerica(expr)(np.linspace(a, b, 2001)))
    return float(np.max(v)) if np.all(np.isfinite(v)) else math.inf


def texto_cota(K, potencia, divisor, a, b, m, n_cifras):
    if not math.isfinite(K):
        return "no hay cota: la derivada que pide la fórmula no está acotada en [a, b]"
    return cifras(K * (b - a)**potencia / (divisor * m**(potencia - 1)), n_cifras)


def tabla_error(f, a, b, exacto, ns=(2, 4, 8, 16, 32)):
    """Filas (n, T_n, error T, razón T, S_n, error S, razón S)."""
    filas, prev = [], None
    for m in ns:
        T, S = reglas_funcion(f, a, b, m)
        eT, eS = T - float(exacto), S - float(exacto)
        rT = abs(prev[0] / eT) if prev and eT else None
        rS = abs(prev[1] / eS) if prev and eS else None
        filas.append((m, T, eT, rT, S, eS, rS)); prev = (eT, eS)
    return filas


def calculadora_trapecio_simpson(modo, f_txt, a, b, m, datos_txt, h, n_cifras):
    try:
        if modo == "función":
            f = parsear(f_txt)
            T, S = reglas_funcion(f, a, b, m)
            K2, K4 = cota_maxima(sp.diff(f, x, 2), a, b), cota_maxima(sp.diff(f, x, 4), a, b)
            print(f"({n_cifras} cifras significativas, redondeado; la tabla de error usa 8)")
            print(f"T_{m} ≈ {cifras(T, n_cifras)}   cota del error: {texto_cota(K2, 3, 12, a, b, m, n_cifras)}")
            if S is None:
                print(f"S_{m}: n impar, Simpson no aplica")
            else:
                print(f"S_{m} ≈ {cifras(S, n_cifras)}   cota del error: {texto_cota(K4, 5, 180, a, b, m, n_cifras)}")
            exacto = integral_exacta(f, a, b)
            if exacto.is_number and not exacto.has(sp.Integral):
                print(f"valor exacto (sympy): {exacto} ≈ {cifras(exacto, n_cifras)}")
                print("   n          T_n      error T   razón        S_n      error S   razón")
                for fila in tabla_error(f, a, b, exacto):
                    mm, T_, eT, rT, S_, eS, rS = fila
                    print(f"{mm:4d}  {cifras(T_, 8):>11}  {eT:10.2e}  {('%.2f' % rT) if rT else '':>6}"
                          f"  {cifras(S_, 8):>11}  {eS:10.2e}  {('%.2f' % rS) if rS else '':>6}")
        else:
            y = lista_numeros(datos_txt, "las lecturas")
            T = trapecio_datos(y, h)            # valida h y el número de lecturas antes de imprimir
            print(f"{len(y)} lecturas, {len(y) - 1} franjas de ancho {h}   ({n_cifras} cifras significativas, redondeado)")
            print(f"trapecio ≈ {cifras(T, n_cifras)}")
            print(f"Simpson ≈ {cifras(simpson_datos(y, h), n_cifras)}")
    except (ValueError, TypeError) as err:
        print("Revisa la entrada:", err)


widgets.interact(calculadora_trapecio_simpson,
    modo=widgets.Dropdown(options=["función", "tabla"], value="función", description="modo"),
    f_txt=widgets.Text(value="1/x", description="f(x) ="),
    a=widgets.FloatText(value=1, description="a"), b=widgets.FloatText(value=2, description="b"),
    m=widgets.IntSlider(value=4, min=1, max=64, description="n"),
    datos_txt=widgets.Text(value="40, 42, 45, 44, 41", description="lecturas"),
    h=widgets.FloatText(value=2, description="ancho h"),
    n_cifras=widgets.IntSlider(value=6, min=1, max=10, description="cifras sig."));

# %% [markdown]
# **Las cuatro reglas lado a lado.** La tabla de la figura 2.8 del libro: error de $L_n$, $M_n$, $T_n$ y $S_n$ para $P(0\le Z\le1)=\int_0^1\varphi(x)\,dx$ y la razón entre errores consecutivos. Al duplicar $n$, las razones se acercan a 2, 4, 4 y 16.

# %%
def tabla_cuatro_reglas(f, a, b, exacto, ns=(2, 4, 8, 16, 32, 64)):
    F, ex = numerica(f), float(exacto)
    filas, previo = [], None
    for m in ns:
        xs = np.linspace(a, b, m + 1); h = (b - a) / m; y = F(xs)
        valores = (h * y[:-1].sum(), h * F((xs[:-1] + xs[1:]) / 2).sum(), trapecio_datos(y, h), simpson_datos(y, h))
        errores = [v - ex for v in valores]
        razones = [abs(p / e) if previo else None for p, e in zip(previo or errores, errores)]
        filas.append((m, errores, razones)); previo = errores
    return filas


_phi = sp.exp(-x**2 / 2) / sp.sqrt(2 * sp.pi)
_tabla4 = tabla_cuatro_reglas(_phi, 0, 1, sp.erf(1 / sp.sqrt(2)) / 2)
print("   n      error L      error M      error T      error S   razones L, M, T, S")
for m, errores, razones in _tabla4:
    r = "" if razones[0] is None else "  ".join(f"{q:5.2f}" for q in razones)
    print(f"{m:4d}  " + "  ".join(f"{e:11.3e}" for e in errores) + "   " + r)

# %%
# Casos de prueba de la sección 5 (resultado conocido)
PRUEBAS_5 = [
    ("x^2 en [0, 1]: S_2 = 1/3 exacto (Simpson es exacto en cúbicas)", lambda: cerca(reglas_funcion(x**2, 0, 1, 2)[1], 1 / 3)),
    ("1/x en [1, 2], n = 4: T ≈ 0.697024, S ≈ 0.693254",
     lambda: cerca(reglas_funcion(1 / x, 1, 2, 4)[0], 0.6970238095, 1e-9) and cerca(reglas_funcion(1 / x, 1, 2, 4)[1], 0.6932539683, 1e-9)),
    ("e^x en [0, 1], n = 2: T ≈ 1.75393, S ≈ 1.71886",
     lambda: cerca(reglas_funcion(sp.exp(x), 0, 1, 2)[0], 1.7539310925, 1e-9) and cerca(reglas_funcion(sp.exp(x), 0, 1, 2)[1], 1.7188611519, 1e-9)),
    ("tabla constante de 6 lecturas: trapecio 25 y Simpson rechazado",
     lambda: cerca(trapecio_datos([5] * 6, 1), 25) and _rechaza(lambda: simpson_datos([5] * 6, 1))),
    ("tabla de las cuatro reglas para la normal: razones ≈ 2, 4, 4 y 16 con n = 64",
     lambda: all(abs(q - r) < t for q, r, t in zip(_tabla4[-1][2], (2, 4, 4, 16), (0.05, 0.05, 0.05, 0.3)))),
    ("razones de la tabla de error para 1/x: ≈ 4 (trapecio) y ≈ 16 (Simpson)",
     lambda: abs(tabla_error(1 / x, 1, 2, sp.log(2))[-1][3] - 4) < 0.05 and abs(tabla_error(1 / x, 1, 2, sp.log(2))[-1][6] - 16) < 0.5),
    ("la cota del trapecio para 1/x con n = 409 es menor que 1e-6",
     lambda: cota_maxima(sp.diff(1 / x, x, 2), 1, 2) / (12 * 409**2) < 1e-6),
    ("a >= b o h <= 0 se rechazan; √x en [0, 1] no tiene cota", lambda: _rechaza(lambda: reglas_funcion(x, 1, 1, 4))
                                                                 and _rechaza(lambda: trapecio_datos([1, 2], 0))
                                                                 and math.isinf(cota_maxima(sp.diff(sp.sqrt(x), x, 2), 0, 1))),
]
for nombre, prueba in PRUEBAS_5:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Retoma tu predicción de la semilla.** La celda siguiente cuenta cuántas franjas necesita cada regla para que $\int_0^1\frac{4}{1+x^2}dx$ dé $\pi$ con seis cifras ($|E|<5\times10^{-6}$). La suma izquierda tarda unos segundos: son muchos rectángulos.

# %%
def _reglas_pi(m, regla):
    xs = np.linspace(0, 1, m + 1); y = 4 / (1 + xs**2); h = 1 / m
    if regla == "L": return h * y[:-1].sum()
    if regla == "M": xm = (xs[:-1] + xs[1:]) / 2; return h * (4 / (1 + xm**2)).sum()
    if regla == "T": return trapecio_datos(y, h)
    return simpson_datos(y, h)


def minimo_n(regla, tol=5e-6):
    paso = 2 if regla == "S" else 1
    lo = paso
    while abs(_reglas_pi(lo, regla) - math.pi) >= tol:     # duplica hasta pasarse
        lo *= 2
    hi, lo = lo, max(paso, lo // 2)
    while hi - lo > paso:                                  # bisección sobre n
        mid = (lo + hi) // 2
        mid -= mid % paso
        if abs(_reglas_pi(mid, regla) - math.pi) < tol: hi = mid
        else: lo = mid
    return hi


for regla, nombre in (("L", "suma izquierda"), ("M", "punto medio"), ("T", "trapecio"), ("S", "Simpson")):
    print(f"{nombre:15s}: n = {minimo_n(regla)}")

# %% [markdown]
# **Contrasta (sección 5).** Calcula a mano $T_2$ y $S_2$ para $\int_0^2x^3\,dx$ y compáralos con el valor exacto, $4$.

# %%
mi_T2 = None    # escribe un número
mi_S2 = None    # escribe un número

if mi_T2 is None or mi_S2 is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    T2, S2 = reglas_funcion(x**3, 0, 2, 2)
    print("La calculadora da: T_2 =", cifras(T2, 6), "  S_2 =", cifras(S2, 6))
    print("coinciden" if cerca(mi_T2, T2, 1e-9) and cerca(mi_S2, S2, 1e-9) else
          "NO coinciden: con Δx = 1, las lecturas son f(0) = 0, f(1) = 1 y f(2) = 8")
