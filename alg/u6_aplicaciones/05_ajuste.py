# ID: ALG-U6-NB05
# Notebook: alg/u6_aplicaciones.ipynb · sección 6.5 ajuste de datos por mínimos cuadrados
# Repositorio: alg/u6_aplicaciones/05_ajuste.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 5. Ajuste de datos por mínimos cuadrados
#
# **Repaso de la Unidad 5.** Con más ecuaciones que incógnitas, $A\boldsymbol\theta=\mathbf y$ casi nunca tiene solución; mínimos cuadrados busca el $\hat{\boldsymbol\theta}$ que hace más pequeño $\|\mathbf y-A\boldsymbol\theta\|$ y lo da la solución de las **ecuaciones normales** $A^TA\,\hat{\boldsymbol\theta}=A^T\mathbf y$. En la práctica se usa `np.linalg.lstsq`, que no forma $A^TA$. Lo único que se pide es que el modelo sea **lineal en los parámetros**: cada dato da un renglón de $A$.
#
# **Recta por el origen.** La ley de Ohm $V=R\,I$, la de Hooke $F=k\,x$ o la dilatación $\Delta L=\alpha L_0\,\Delta T$ tienen la forma $y=k\,x$: $A$ es una sola columna (los $x_i$) y las ecuaciones normales se reducen a
#
# $$k=\frac{\sum_i x_iy_i}{\sum_i x_i^2}.$$
#
# **Modelo discreto de primer orden.** Un horno, un motor o un tanque que se mide cada $\Delta t$ segundos se describe muchas veces con
#
# $$x[k+1]=a\,x[k]+b\,u[k],$$
#
# donde $u[k]$ es la entrada (el voltaje, la apertura de una válvula) y $x[k]$ es la **desviación de la salida respecto al reposo** (por ejemplo, la temperatura del horno sobre la del ambiente, no la temperatura absoluta). Como el modelo no tiene término constante, con la temperatura absoluta el ajuste queda sesgado: el modelo diría que sin entrada la salida tiende a 0 °C. Es lineal en $a$ y $b$: **cada paso medido es un renglón** $\big(x[k],\ u[k]\big)$ de $A$, con $x[k+1]$ del lado derecho. Con $N$ entradas hacen falta $N+1$ salidas, $x[0],\dots,x[N]$. De $a$ y $b$ salen dos constantes físicas, solo si $0<a<1$ (el sistema se estabiliza sin oscilar):
#
# $$\tau=-\frac{\Delta t}{\ln a}\ \ (\text{constante de tiempo}),\qquad K=\frac{b}{1-a}\ \ (\text{cuánto sube la salida en estado estacionario por unidad de entrada}).$$
#
# Si $a\le 0$ o $a\ge 1$ no hay constante de tiempo: con $a=0$ la salida copia la entrada un paso después ($K=b$); con $-1<a<0$ converge a $K\,u$ oscilando, pero $\tau$ no existe porque $a=e^{-\Delta t/\tau}$ tendría que ser positivo; con $a\le-1$ o $a\ge 1$ no se estabiliza.
#
# **Validar.** La predicción a un paso usa el dato medido $x[k]$ para predecir $x[k+1]$, así que casi siempre se ve bien. La **simulación libre** es una prueba más exigente: arrancar en $x[0]$ y dejar que el modelo use sus propias predicciones; si se aleja de los datos, el modelo no sirve aunque el error a un paso sea chico. Mejor todavía es comparar el modelo con datos que no se usaron en el ajuste (por ejemplo, otra corrida del mismo equipo).
#
# **Columnas casi dependientes.** Si una columna de $A$ es casi un múltiplo (o una combinación) de otras, el número de condición de $A$ es grande: el ajuste puede verse bien, pero los coeficientes **por separado** no son confiables (al quitar un dato cambian mucho), y las predicciones solo valen cerca de los datos. La calculadora avisa cuando el número de condición pasa de 1000.

# %%
def _pares(x, y, nx, ny):
    xs, ys = datos(x, nx), datos(y, ny)
    if len(xs) != len(ys):
        raise ValueError(f"hay {_cuenta(len(xs), 'dato', 'datos')} de x y {len(ys)} de y: debe haber un y por cada x"
                         f"{pista_coma(x) or pista_coma(y) or _pista_decimal(x) or _pista_decimal(y)}")
    return xs, ys


def ajuste_por_origen(x, y):
    """Recta por el origen y = k x por mínimos cuadrados: k = Σxy/Σx². Exacto con hasta 12 datos; con más, en
    decimales. Diccionario con k, residuos (y - kx), SSE, RMSE, Sxy, Sxx, exacto, x, y y n."""
    xs, ys = _pares(x, y, "los datos de x", "los datos de y")
    if all(v == 0 for v in xs):
        raise ValueError("todos los x son 0: la recta por el origen no queda determinada")
    n = len(xs)
    exacto = n <= _MAX_EXACTO
    if exacto:
        sxy, sxx = sum(a * b for a, b in zip(xs, ys)), sum(a * a for a in xs)
        k = sxy / sxx
        res = [b - k * a for a, b in zip(xs, ys)]
        sse = sum(r ** 2 for r in res)
        rmse = sp.sqrt(sse / n)
    else:
        X, Y = np.array([float(v) for v in xs]), np.array([float(v) for v in ys])
        sxy, sxx = float(X @ Y), float(X @ X)
        k = sxy / sxx
        res = list(Y - k * X)
        sse = float(np.dot(res, res))
        rmse = math.sqrt(sse / n)
    return {"k": k, "residuos": res, "SSE": sse, "RMSE": rmse, "Sxy": sxy, "Sxx": sxx, "exacto": exacto,
            "x": xs, "y": ys, "n": n}


def ajuste_lineal(columnas, y, nombres=None):
    """Ajuste por mínimos cuadrados de y ≈ θ1·c1 + θ2·c2 + ... con las columnas dadas (listas de números del mismo
    largo que y), en decimales con np.linalg.lstsq. Diccionario con coeficientes, nombres, residuos, RMSE,
    condicion (número de condición de A) y avisos: si condicion pasa de 1000, «columnas casi dependientes: los
    coeficientes no son confiables por separado». Si las columnas son dependientes, ValueError."""
    if isinstance(columnas, (str, np.ndarray)) or not isinstance(columnas, (list, tuple)) or not columnas:
        raise ValueError("escribe las columnas como una lista de listas de números, una lista por columna")
    if len(columnas) > 6:
        raise ValueError(f"escribiste {len(columnas)} columnas; esta función acepta hasta 6")
    cols = [datos(c, f"la columna {j}") for j, c in enumerate(columnas, 1)]
    ys = datos(y, "los datos de y")
    nombres = list(nombres) if nombres is not None else [f"θ{j}" for j in range(1, len(cols) + 1)]
    if len(nombres) != len(cols):
        raise ValueError(f"hay {_cuenta(len(cols), 'columna', 'columnas')} y {_cuenta(len(nombres), 'nombre', 'nombres')}: "
                         "debe haber un nombre por columna")
    for j, c in enumerate(cols, 1):
        if len(c) != len(ys):
            raise ValueError(f"la columna {j} tiene {_cuenta(len(c), 'dato', 'datos')} y y tiene {len(ys)}: deben ser "
                             "los mismos")
    if len(ys) < len(cols):
        raise ValueError(f"hay {_cuenta(len(ys), 'dato', 'datos')} y {len(cols)} coeficientes: hacen falta al menos "
                         "tantos datos como coeficientes")
    A = np.column_stack([[float(v) for v in c] for c in cols])
    Y = np.array([float(v) for v in ys])
    with np.errstate(all="ignore"):
        rango = int(np.linalg.matrix_rank(A))
        if rango < A.shape[1]:
            raise Diagnostico(f"las columnas son dependientes (rango {rango} < {A.shape[1]}): los coeficientes no quedan "
                             "determinados")
        cond = float(np.linalg.cond(A))
        coef = np.linalg.lstsq(A, Y, rcond=None)[0]
    pred = A @ coef
    res = Y - pred
    avisos = []
    if cond > _COND_AVISO:
        avisos.append("columnas casi dependientes: los coeficientes no son confiables por separado")
    return {"coeficientes": [float(c) for c in coef], "nombres": nombres, "residuos": [float(r) for r in res],
            "predichos": [float(p) for p in pred], "RMSE": math.sqrt(float(res @ res) / len(Y)), "condicion": cond,
            "avisos": avisos, "A": A}


def _decimal(v, nombre):
    """Un número como decimal de Python: el texto se lee con numero(); un número (de Python, numpy o sympy) se
    convierte y debe ser finito. Mensajes en español."""
    if isinstance(v, str):
        v = numero(v, nombre)
    try:
        f = float(v)
    except (TypeError, ValueError):
        raise ValueError(f"{nombre} debe ser un número ({_AYUDA_NUMERO})") from None
    if not math.isfinite(f):
        raise ValueError(f"{nombre} debe ser un número finito (escribiste {v})")
    return f


def simula_arx(a, b, u, x0):
    """Simulación libre de x[k+1] = a·x[k] + b·u[k] desde x[0] = x0: devuelve [x[0], x[1], ..., x[N]] (N = len(u)).
    Exacta (fracciones) con hasta 12 pasos y datos racionales; con más, en decimales. Si la simulación diverge, los
    valores que ya no caben en los decimales quedan como inf (identifica_arx lo avisa)."""
    us = datos(u, "la entrada u")
    if len(us) <= _MAX_EXACTO:
        a, b, x = numero(a, "a"), numero(b, "b"), numero(x0, "x0")
        salida = [x]
        for uk in us:
            x = a * x + b * uk
            salida.append(x)
        return salida
    a, b, x = _decimal(a, "a"), _decimal(b, "b"), _decimal(x0, "x0")
    salida = [x]
    for uk in us:
        x = a * x + b * float(uk) if math.isfinite(x) else math.inf
        salida.append(x if math.isfinite(x) else math.inf)
    return salida


def identifica_arx(u, x, dt, n_cifras=4):
    """Identifica a y b del modelo x[k+1] = a·x[k] + b·u[k] por mínimos cuadrados: cada paso es un renglón (x[k], u[k])
    de A y x[k+1] va del lado derecho. Necesita len(x) = len(u) + 1. Exacto con hasta 12 pasos; con más, en decimales
    con np.linalg.lstsq. Diccionario con a, b, tau = -Δt/ln a y K = b/(1 - a) (solo si 0 < a < 1; si no, None y un
    aviso; con a = 0, K = b y con -1 < a < 0, K = b/(1 - a), sin τ), A, y, prediccion (a un paso), residuos, RMSE_1 (a un
    paso), simulacion (libre, desde x[0]), RMSE_libre, cond (número de condición de A), avisos, exacto, u, x, dt y n.
    n_cifras: las cifras de los decimales en los avisos."""
    us, xs = datos(u, "la entrada u"), datos(x, "la salida x")
    if len(xs) != len(us) + 1:
        raise ValueError(f"x debe tener un dato más que u (x[0], ..., x[N] y u[0], ..., u[N-1]): escribiste "
                         f"{_cuenta(len(us), 'dato', 'datos')} de u y {len(xs)} de x"
                         f"{pista_coma(u) or pista_coma(x) or _pista_decimal(u) or _pista_decimal(x)}")
    n = len(us)
    if n < 2:
        raise ValueError("hacen falta al menos dos pasos (dos ecuaciones) para encontrar a y b: escribe al menos 2 datos "
                         "de u y 3 de x")
    h = numero(dt, "Δt")
    if h <= 0:
        raise ValueError(f"Δt debe ser positivo (escribiste {_ntxt(h)})")
    exacto = n <= _MAX_EXACTO
    dependientes = ("las columnas x[k] y u[k] de A son dependientes (por ejemplo, x es proporcional a u en todos los "
                    "pasos): los datos no distinguen a de b; usa una entrada que cambie (escalones) o más pasos")
    if exacto:
        A, Y = sp.Matrix([[xs[k], us[k]] for k in range(n)]), sp.Matrix(xs[1:])
        G = A.T * A
        if G.det() == 0:
            raise Diagnostico(dependientes)
        a, b = G.LUsolve(A.T * Y)
        pred = list(A * sp.Matrix([a, b]))
        res = [v - p for v, p in zip(xs[1:], pred)]
        rmse1 = sp.sqrt(sum(r ** 2 for r in res) / n)
        cond = float(np.linalg.cond(_flotante(A)))
    else:
        A = np.column_stack([[float(v) for v in xs[:-1]], [float(v) for v in us]])
        Y = np.array([float(v) for v in xs[1:]])
        with np.errstate(all="ignore"):
            if int(np.linalg.matrix_rank(A)) < 2:
                raise Diagnostico(dependientes)
            cond = float(np.linalg.cond(A))
            a, b = (float(c) for c in np.linalg.lstsq(A, Y, rcond=None)[0])
        pred = list(A @ np.array([a, b]))
        res = [float(v) - p for v, p in zip(Y, pred)]
        rmse1 = math.sqrt(sum(r * r for r in res) / n)
    sim = simula_arx(a, b, us, xs[0])
    diverge = not all(math.isfinite(float(s)) and abs(float(s)) < 1e100 for s in sim)
    dif = [v - s for v, s in zip(xs[1:], sim[1:])]
    if diverge:
        rmse_libre = math.inf
    elif exacto:
        rmse_libre = sp.sqrt(sum(d ** 2 for d in dif) / n)
    else:
        rmse_libre = math.sqrt(sum(float(d) ** 2 for d in dif) / n)
    avisos, tau, K = [], None, None
    if 0 < a < 1:
        K = b / (1 - a)
        tau = -h / sp.log(a) if exacto else -float(h) / math.log(a)
    elif a == 0:
        K = b
        avisos.append("a = 0: la salida copia la entrada un paso después (x[k+1] = b·u[k]), sin dinámica; K = b y no "
                      "hay constante de tiempo τ")
    elif -1 < a < 0:
        K = b / (1 - a)
        avisos.append("-1 < a < 0: la salida converge a K·u oscilando (cambia de lado en cada paso); K = b/(1 - a) "
                      "existe, pero no hay constante de tiempo τ, porque a = e^(-Δt/τ) tendría que ser positivo")
    elif a <= -1:
        avisos.append("a ≤ -1: la salida oscila sin amortiguarse o crece en cada paso; el modelo no se estabiliza, así "
                      "que no hay constante de tiempo τ ni ganancia K")
    else:
        avisos.append("a ≥ 1: el modelo no se estabiliza (la salida crece sin límite o no regresa), así que no hay "
                      "constante de tiempo τ ni ganancia K")
    if cond > _COND_AVISO:
        avisos.append(f"el número de condición de A es ≈ {cifras(cond, n_cifras)}: columnas casi dependientes, a y b "
                      "no son confiables por separado")
    escala = max([abs(float(v)) for v in xs] + [1e-300])
    if diverge:
        avisos.append("la simulación libre diverge (sus valores crecen sin límite): el modelo no reproduce el sistema")
    elif float(rmse_libre) > 5 * float(rmse1) and float(rmse_libre) > 1e-6 * escala:
        avisos.append("la simulación libre se aleja de los datos mucho más que la predicción a un paso: el modelo no "
                      "reproduce bien el sistema")
    return {"a": a, "b": b, "tau": tau, "K": K, "A": A, "y": Y, "prediccion": pred, "residuos": res, "RMSE_1": rmse1,
            "simulacion": sim, "RMSE_libre": rmse_libre, "cond": cond, "avisos": avisos, "exacto": exacto, "u": us,
            "x": xs, "dt": h, "n": n}


def _informe_origen(res, k_nombre, unidad_k, n_cifras):
    print(f"{_cuenta(res['n'], 'dato', 'datos')}; modelo y = {k_nombre}·x (A es la columna de los x).")
    if res["exacto"]:
        print(f"Σ x·y = {_ntxt(res['Sxy'])},  Σ x² = {_ntxt(res['Sxx'])}  →  {k_nombre} = Σxy/Σx² = "
              f"{_cantidad(res['k'], unidad_k, n_cifras)}")
        print(f"Residuos y - {k_nombre}·x = {_tupla_u(sp.Matrix(res['residuos']), '', n_cifras)}")
        print(f"{_magnitud('SSE', res['SSE'], '', n_cifras)},  {_magnitud('RMSE', res['RMSE'], '', n_cifras)}")
    else:
        print(f"{k_nombre} ≈ {cifras(res['k'], n_cifras)} {unidad_k};  RMSE ≈ {cifras(res['RMSE'], n_cifras)}")


def _informe_arx(res, ut, n_cifras):
    n, a, b = res["n"], res["a"], res["b"]
    nt = lambda v: _ntxt(v, cif=n_cifras)
    print(f"Modelo x[k+1] = a·x[k] + b·u[k];  N = {n} pasos (una ecuación por paso), Δt = {nt(res['dt'])} {ut}.")
    if n <= 6:
        for k in range(n):
            print(f"  k = {k}:  x[{k + 1}] = a·x[{k}] + b·u[{k}]  →  {nt(res['x'][k + 1])} = "
                  f"{nt(res['x'][k])}·a + {nt(res['u'][k])}·b")
    if res["exacto"]:
        print("Cálculo exacto con fracciones: ecuaciones normales AᵀA θ = Aᵀy, con θ = (a, b).")
        if n <= 6:
            muestra("A (columnas x[k] y u[k])", res["A"], n_cifras, tex=r"A\ (\text{columnas } x[k] \text{ y } u[k])")
    else:
        print(f"Cálculo en decimales con np.linalg.lstsq, porque hay más de {_MAX_EXACTO} pasos.")
    print(f"{_magnitud('a', a, '', n_cifras)},  {_magnitud('b', b, '', n_cifras)}")
    print(f"Número de condición de A ≈ {cifras(res['cond'], n_cifras)}")
    if res["K"] is not None:
        print(_magnitud("K = b/(1 - a)" if a != 0 else "K = b", res["K"], "", n_cifras) +
              "  (subida de la salida en estado estacionario por unidad de entrada)")
    if res["tau"] is not None:
        print(f"τ = -Δt/ln a = -{nt(res['dt'])}/ln({_dec_o_cifras(a, n_cifras)}) ≈ {cifras(res['tau'], n_cifras)} {ut}  "
              "(constante de tiempo)")
    medida = lambda v: "= 0" if (isinstance(v, sp.Basic) and v == 0) else f"≈ {cifras(v, n_cifras)}"
    print(f"RMSE a un paso {medida(res['RMSE_1'])};  RMSE de la simulación libre {medida(res['RMSE_libre'])}" +
          ("  (el modelo pasa exactamente por todos los datos)" if res["RMSE_libre"] == 0 else ""))
    for aviso in res["avisos"]:
        print("Aviso:", aviso + ".")


def dibuja_identificacion(resultado, unidad_t="s", n_cifras=4):
    """Tres paneles: la salida medida con la simulación libre y la predicción a un paso; la entrada en escalones; y
    los residuos a un paso."""
    n, h = resultado["n"], float(resultado["dt"])
    t = np.arange(n + 1) * h
    X = np.array([float(v) for v in resultado["x"]])
    S = np.array([float(v) for v in resultado["simulacion"]])
    S[~np.isfinite(S) | (np.abs(S) > 1e100)] = np.nan             # una simulación que diverge no se dibuja completa
    P = np.array([float(v) for v in resultado["prediccion"]])
    U = np.array([float(v) for v in resultado["u"]])
    R = np.array([float(v) for v in resultado["residuos"]])
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8.6, 8.4), sharex=True,
                                        gridspec_kw={"height_ratios": [2.3, 1.0, 1.2]})
    marca = n <= 60
    ax1.plot(t, X, "o", color="black", ms=6 if marca else 3, label="x medida", zorder=3)
    rmse = lambda v: "= 0" if (isinstance(v, sp.Basic) and v == 0) else f"≈ {cifras(v, n_cifras)}"
    ax1.plot(t, S, "--", color=_COLORES[0], lw=2.2, marker="s" if marca else None, ms=4,
             label=f"simulación libre desde x[0] (RMSE {rmse(resultado['RMSE_libre'])})")
    ax1.plot(t[1:], P, ":", color=_COLORES[1], lw=2.0, marker="^" if marca else None, ms=5,
             label=f"predicción a un paso (RMSE {rmse(resultado['RMSE_1'])})")
    ax1.set_ylabel("salida x"); ax1.grid(alpha=0.4); ax1.legend(fontsize=8, loc="best")
    igual = lambda v: ("= " if isinstance(v, sp.Basic) and _decimal_exacto(v, n_cifras) else "≈ ") + _dec_o_cifras(v, n_cifras)
    ax1.set_title(f"x[k+1] = a·x[k] + b·u[k] con a {igual(resultado['a'])} y b {igual(resultado['b'])}"
                  + (f";  τ ≈ {cifras(resultado['tau'], n_cifras)} {unidad_t}" if resultado["tau"] is not None else ""),
                  fontsize=10)
    ax2.step(t, np.append(U, U[-1]), where="post", color=_COLORES[2], lw=2.0, label="entrada u (escalones)")
    ax2.set_ylabel("entrada u"); ax2.grid(alpha=0.4); ax2.legend(fontsize=8, loc="best")
    escala = max(float(np.max(np.abs(X))), 1e-300)
    if np.all(np.abs(R) <= 1e-12 * escala):
        R = np.zeros_like(R)
    ax3.vlines(t[1:], 0, R, colors="0.35", lw=1.6)
    ax3.plot(t[1:], R, "s", color=_COLORES[3], ms=5, label="residuo a un paso x[k+1] - predicción")
    ax3.axhline(0, color="black", lw=1.0, ls="--")
    ax3.set_ylabel("residuo"); ax3.set_xlabel(f"tiempo ({unidad_t})"); ax3.grid(alpha=0.4)
    ax3.legend(fontsize=8, loc="best")
    plt.tight_layout()
    plt.show()

# %% [markdown]
# **Recorrido paso a paso.** (1) Ley de Ohm por el origen: con $I=0.1,\ 0.2,\ 0.3,\ 0.4$ A se miden $V=1.9,\ 4.1,\ 6.0,\ 8.1$ V; $\sum IV=6.05$ y $\sum I^2=0.30$, así que $R=6.05/0.30=\tfrac{121}{6}\approx 20.17\ \Omega$. (2) Modelo de primer orden con $u=1,1,1$ y $x=0,\ 2,\ 3.6,\ 4.88$: la primera ecuación, $2=a\cdot0+b\cdot1$, da $b=2$, y la segunda, $3.6=2a+2$, da $a=0.8$; la tercera también se cumple, así que $a$ y $b$ son exactos. Entonces $K=2/(1-0.8)=10$ y, con $\Delta t=1$ s, $\tau=-1/\ln 0.8\approx 4.48$ s. (3) Una serie de 12 pasos generada con $a=0.95$, $b=0.5$, la entrada $u=1,1,1,0,0,0,0,0,0,1,1,1$ y $x[0]=0$: la identificación recupera $a$ y $b$, y la simulación libre pasa por todos los datos. (4) Columnas casi dependientes: con $x_2\approx 10\,x_1$, el ajuste se ve bien, pero al quitar un solo dato los coeficientes cambian mucho.

# %%
print("(1) Ley de Ohm, V = R·I")
_ohm = ajuste_por_origen("0.1 0.2 0.3 0.4", "1.9 4.1 6.0 8.1")
_informe_origen(_ohm, "R", "Ω", 4)
print("\n" + "=" * 90 + "\n\n(2) Modelo de primer orden con datos exactos")
_arx = identifica_arx("1 1 1", "0 2 3.6 4.88", 1)
_informe_arx(_arx, "s", 4)
dibuja_identificacion(_arx, "s")
print("=" * 90 + "\n\n(3) Serie sintética de 12 pasos con a = 0.95 y b = 0.5")
_u_sint = [1, 1, 1, 0, 0, 0, 0, 0, 0, 1, 1, 1]
_x_sint = simula_arx(0.95, 0.5, _u_sint, 0)
print("x generada:", ", ".join(cifras(v, 4) for v in _x_sint))
_sint = identifica_arx(_u_sint, _x_sint, 1)
_informe_arx(_sint, "s", 4)
dibuja_identificacion(_sint, "s")

# %%
print("(4) Columnas casi dependientes: x2 = 10·x1 + (±0.01)")
_x1 = [1, 2, 3, 4, 5, 6]
_x2 = [10 * a + d for a, d in zip(_x1, [0.01, -0.01, 0.01, -0.01, 0.01, -0.01])]
_yc = [3.0, 6.2, 8.9, 12.1, 15.1, 17.6]
_todos = ajuste_lineal([_x1, _x2], _yc, ["θ1", "θ2"])
_sin_ultimo = ajuste_lineal([_x1[:-1], _x2[:-1]], _yc[:-1], ["θ1", "θ2"])
for _nombre, _r in (("con los 6 datos", _todos), ("sin el último dato", _sin_ultimo)):
    print(f"  {_nombre}: θ1 ≈ {cifras(_r['coeficientes'][0], 4)}, θ2 ≈ {cifras(_r['coeficientes'][1], 4)};  "
          f"RMSE ≈ {cifras(_r['RMSE'], 4)};  número de condición ≈ {cifras(_r['condicion'], 4)}")
print("  Aviso:", _todos["avisos"][0] + ".")
print(f"  La combinación θ1 + 10·θ2, que es lo que los datos sí determinan, casi no cambia: "
      f"{cifras(_todos['coeficientes'][0] + 10 * _todos['coeficientes'][1], 4)} y "
      f"{cifras(_sin_ultimo['coeficientes'][0] + 10 * _sin_ultimo['coeficientes'][1], 4)}.")

# %% [markdown]
# **Calculadora.** Escribe la entrada $u[0],\dots,u[N-1]$ y la salida $x[0],\dots,x[N]$ (un dato más que $u$), separadas por espacios o comas, el paso de tiempo $\Delta t$ y su unidad. La salida es la **desviación respecto al reposo** (por ejemplo, la temperatura sobre la del ambiente): como el modelo no tiene término constante, con valores absolutos el ajuste queda sesgado. Con hasta 12 pasos el cálculo es exacto; con más (hasta 500), en decimales. Prueba agregarle ruido a la serie sintética (cambia algunos datos en la segunda cifra decimal), una entrada constante con la salida también constante (columnas dependientes) y una salida que crece sin parar (por ejemplo `0 1 2.2 3.64` con `1 1 1`).

# %%
def calculadora_arx(u_txt, x_txt, dt_txt, unidad_t, n_cifras):
    with bloque():                 # el texto aparece completo, de una sola vez; las gráficas van después
        try:
            ut = nombre_unidad(unidad_t, "tiempo", "s o min")
            res = identifica_arx(u_txt, x_txt, dt_txt, n_cifras=n_cifras)
        except ValueError as err:
            _reporta(err); return
        _aviso_coma(u_txt, x_txt)
        _informe_arx(res, ut, n_cifras)
    dibuja_identificacion(res, ut, n_cifras)


_CONTROLES_5 = dict(
    u_txt=widgets.Textarea(value="1 1 1", description="u (entrada):", continuous_update=False,
                           layout=widgets.Layout(width="420px", height="60px"), style={"description_width": "initial"}),
    x_txt=widgets.Textarea(value="0 2 3.6 4.88", description="x (salida):", continuous_update=False,
                           layout=widgets.Layout(width="420px", height="60px"), style={"description_width": "initial"}),
    dt_txt=widgets.Text(value="1", description="Δt =", continuous_update=False),
    unidad_t=widgets.Text(value="s", description="unidad de tiempo:", continuous_update=False,
                          style={"description_width": "initial"}),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_arx, **_CONTROLES_5);

# %%
# Casos de prueba de la sección 5 (resultado conocido)
PRUEBAS_5 = [
    ("ley de Ohm por el origen con I = 0.1, ..., 0.4 A y V = 1.9, 4.1, 6.0, 8.1 V → R = 121/6 (tolerancia 1e-12)",
     lambda: _cerca(ajuste_por_origen("0.1 0.2 0.3 0.4", "1.9 4.1 6.0 8.1")["k"], 121 / 6, 1e-12)
             and ajuste_por_origen("0.1 0.2 0.3 0.4", "1.9 4.1 6.0 8.1")["k"] == sp.Rational(121, 6)),
    ("ARX exacto con u = 1, 1, 1 y x = 0, 2, 3.6, 4.88 → a = 0.8, b = 2, K = 10, τ = -1/ln 0.8",
     lambda: (lambda r: r["a"] == sp.Rational(4, 5) and r["b"] == 2 and r["K"] == 10
                        and _cerca(r["tau"], -1 / math.log(0.8), 1e-12) and r["RMSE_1"] == 0)(
         identifica_arx("1 1 1", "0 2 3.6 4.88", 1))),
    ("serie sintética (a = 0.95, b = 0.5, 12 pasos) → a y b recuperados (1e-9) y simulación libre con RMSE < 1e-9; "
     "también con 24 pasos (u repetida dos veces), que usa la rama en decimales de simula_arx, identifica_arx y "
     "ajuste_por_origen",
     lambda: (lambda r: _cerca(r["a"], 0.95) and _cerca(r["b"], 0.5) and float(r["RMSE_libre"]) < 1e-9)(
         identifica_arx(_u_sint, simula_arx(0.95, 0.5, _u_sint, 0), 1))
             and (lambda x24: all(isinstance(v, float) for v in x24)
                              and all(_cerca(p, q, 1e-12) for p, q in zip(x24, _x_sint))
                              and (lambda r: not r["exacto"] and _cerca(r["a"], 0.95) and _cerca(r["b"], 0.5)
                                             and r["RMSE_libre"] < 1e-9 and _cerca(r["K"], 10, 1e-8))(
                                  identifica_arx(_u_sint * 2, x24, 1)))(simula_arx(0.95, 0.5, _u_sint * 2, 0))
             and (lambda r: not r["exacto"] and _cerca(r["k"], 2.5) and r["RMSE"] < 1e-9)(
                 ajuste_por_origen(list(range(1, 25)), [2.5 * v for v in range(1, 25)]))),
    ("simula_arx(0.8, 2, [1, 1, 1], 0) = [0, 2, 3.6, 4.88]",
     lambda: simula_arx(0.8, 2, [1, 1, 1], 0) == datos("0 2 3.6 4.88")),
    ("ajuste_lineal con x2 = 10·x1 ± 0.01 → aviso de condición; con columnas independientes, sin aviso",
     lambda: _todos["condicion"] > 1000 and _todos["avisos"] == ["columnas casi dependientes: los coeficientes no son "
                                                                 "confiables por separado"]
             and (lambda r: r["condicion"] < 1000 and not r["avisos"])(
                 ajuste_lineal([_x1, [3, 1, 4, 1, 5, 9]], _yc))),
    ("identifica_arx con len(x) ≠ len(u) + 1 → ValueError «x debe tener un dato más que u» (y con «3,6» la pista del "
     "punto decimal)",
     lambda: _mensaje(lambda: identifica_arx("1 1 1", "0 2 3.6", 1)).startswith("x debe tener un dato más que u")
             and _mensaje(lambda: identifica_arx("1 1", "0 2 3.6 4.88", 1)).startswith("x debe tener un dato más que u")
             and "punto decimal" in _mensaje(lambda: identifica_arx("1 1 1", "0 2 3,6 4,88", 1))),
    ("a ≥ 1 (u = 1, 1, 1 y x = 0, 1, 2.2, 3.64, con a = 1.2) → aviso de que no se estabiliza, sin τ ni K; a = 0 → K = b "
     "sin τ; a = -1/2 → K = 2/3 sin τ (converge oscilando); a = -1 → no se estabiliza",
     lambda: (lambda r: r["a"] == sp.Rational(6, 5) and r["tau"] is None and r["K"] is None
                        and any("no se estabiliza" in a for a in r["avisos"]))(identifica_arx("1 1 1", "0 1 2.2 3.64", 1))
             and (lambda r: r["a"] == 0 and r["K"] == r["b"] == 2 and r["tau"] is None)(
                 identifica_arx("1 1 1", "0 2 2 2", 1))
             and (lambda r: r["a"] == sp.Rational(-1, 2) and r["K"] == sp.Rational(2, 3) and r["tau"] is None
                            and any("oscilando" in a for a in r["avisos"]))(identifica_arx("1 1 1", "0 1 0.5 0.75", 1))
             and (lambda r: r["a"] == -1 and r["K"] is None and r["tau"] is None
                            and any("no se estabiliza" in a for a in r["avisos"]))(identifica_arx("1 1 1", "0 1 0 1", 1))),
    ("prueba de humo: la calculadora corre con sus valores iniciales y con la serie sintética, sin error y sin "
     "«Revisa la entrada»",
     lambda: _humo(calculadora_arx, _CONTROLES_5)
             and _humo(calculadora_arx, _CONTROLES_5, u_txt=" ".join(str(v) for v in _u_sint),
                       x_txt=" ".join(str(v) for v in _x_sint), dt_txt="0.5", unidad_t="min", n_cifras=6)),
]
corre_pruebas(PRUEBAS_5, 5)

# %% [markdown]
# **Contrasta (sección 5).** Con $u=1,\ 1$ y $x=0,\ 3,\ 4.5$, plantea las dos ecuaciones de $x[k+1]=a\,x[k]+b\,u[k]$ (una por paso), resuélvelas a mano y escribe `mis_parametros = (a, b)`.

# %%
mis_parametros = None     # dos números (a, b), por ejemplo (0.4, 2.5)

if mis_parametros is None:     # la referencia se muestra después de tu cálculo
    print("Falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda.")
else:
    try:
        mios = vector(mis_parametros, "mis_parametros")
    except ValueError as err:
        print(f"Revisa tu respuesta: {str(err).rstrip('.')}.")
    else:
        if mios.rows != 2:          # respuesta imposible: la referencia todavía no se muestra
            print(f"El modelo tiene dos parámetros, a y b (escribiste {mios.rows}): revisa tu respuesta y vuelve a "
                  "ejecutar la celda.")
        else:
            ref = identifica_arx("1 1", "0 3 4.5", 1)          # los datos del enunciado
            print(f"La calculadora da: {_magnitud('a', ref['a'], '', 4)}, {_magnitud('b', ref['b'], '', 4)}.")
            if abs(float(mios[0] - ref["a"])) <= 0.01 and abs(float(mios[1] - ref["b"])) <= 0.01:
                print("Coinciden (con tolerancia de 0.01).")
            elif abs(float(mios[1] - ref["a"])) <= 0.01 and abs(float(mios[0] - ref["b"])) <= 0.01:
                print("NO coinciden: probablemente invertiste el orden; a multiplica a x[k] y b a u[k].")
            else:
                print("NO coinciden: el paso k = 0 da x[1] = a·x[0] + b·u[0], con x[0] = 0, así que sale b directo; "
                      "después, el paso k = 1 da a.")
            print("\nEl procedimiento completo:")
            _informe_arx(ref, "s", 4)
