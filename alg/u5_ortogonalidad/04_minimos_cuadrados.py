# ID: ALG-U5-NB04
# Notebook: alg/u5_ortogonalidad.ipynb · sección 5.4 mínimos cuadrados y ajuste de modelos lineales
# Repositorio: alg/u5_ortogonalidad/04_minimos_cuadrados.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 4. Mínimos cuadrados y ajuste de modelos lineales
#
# **El problema.** Con más ecuaciones que incógnitas, un sistema $A\mathbf x=\mathbf b$ casi nunca tiene solución: $\mathbf b$ no está en $\operatorname{Col}(A)$. *Mínimos cuadrados* busca el $\hat{\mathbf x}$ que hace más pequeño el error $\|\mathbf b-A\mathbf x\|$. Por la sección 2, el punto de $\operatorname{Col}(A)$ más cercano a $\mathbf b$ es su proyección: $A\hat{\mathbf x}=\operatorname{proy}_{\operatorname{Col}A}\mathbf b$, y el residuo $\mathbf r=\mathbf b-A\hat{\mathbf x}$ es perpendicular a cada columna de $A$, es decir, $A^T(\mathbf b-A\hat{\mathbf x})=\mathbf 0$. Esas son las **ecuaciones normales**
#
# $$A^TA\,\hat{\mathbf x}=A^T\mathbf b ,$$
#
# con solución única si las columnas de $A$ son independientes.
#
# **Ajuste de una recta.** Con datos $(x_i,y_i)$, $i=1,\dots,n$, el modelo $y=\beta_0+\beta_1x$ da un renglón $\beta_0+\beta_1x_i=y_i$ por dato: $A$ tiene una columna de unos y la columna de los $x_i$, y las ecuaciones normales son $A^TA\,\boldsymbol\beta=A^T\mathbf y$. Se resuelven en forma cerrada:
#
# $$\beta_1=\frac{S_{xy}}{S_{xx}},\qquad \beta_0=\bar y-\beta_1\bar x,\qquad\text{con}\quad S_{xy}=\sum_i(x_i-\bar x)(y_i-\bar y),\quad S_{xx}=\sum_i(x_i-\bar x)^2 .$$
#
# La recta pasa por $(\bar x,\bar y)$: el renglón de la columna de unos en $A^T\mathbf r=\mathbf 0$ dice que **los residuos suman cero**.
#
# **Modelos lineales en los parámetros.** Lo que debe ser lineal son los coeficientes, no la gráfica: un polinomio $y=\beta_0+\beta_1x+\dots+\beta_gx^g$ (columnas $1,x,\dots,x^g$), una señal $y=\beta_0+\beta_1\operatorname{sen}(\omega t)+\beta_2\cos(\omega t)$ con $\omega$ conocida (columnas $1$, $\operatorname{sen}\omega t$, $\cos\omega t$) o un plano $z=\beta_0+\beta_1x+\beta_2y$ (columnas $1,x,y$).
#
# **Linealización.** $y=a\,e^{kx}$ se vuelve $\ln y=\ln a+kx$, una recta en $\ln y$ (necesita $y>0$); $y=a\,x^b$ se vuelve $\ln y=\ln a+b\ln x$ (necesita $x>0$ y $y>0$). Así se minimiza el error en $\ln y$, **no** en $y$: es otro problema de mínimos cuadrados, y los datos con $y$ pequeño pesan más.
#
# **Medir el ajuste.** Con el residuo $\mathbf r=\mathbf y-\hat{\mathbf y}$ (dato menos predicho):
#
# $$\mathrm{SSE}=\|\mathbf r\|^2,\qquad \mathrm{RMSE}=\sqrt{\mathrm{SSE}/n},\qquad R^2=1-\frac{\mathrm{SSE}}{\mathrm{SST}},\quad \mathrm{SST}=\sum_i(y_i-\bar y)^2 .$$
#
# RMSE está en las unidades de $y$. $R^2$ solo tiene sentido si el modelo tiene ordenada al origen (columna de unos); sin ella puede salir negativo, y la calculadora no lo da.
#
# **Precauciones.**
#
# - Un polinomio de grado $n-1$ pasa por los $n$ datos (residuo cero) y no por eso predice mejor: entre los datos y fuera de ellos puede oscilar mucho.
# - Extrapolar (predecir fuera del intervalo de los datos) es riesgoso: nada garantiza que el modelo siga valiendo ahí.
# - Un dato atípico pesa mucho, porque su residuo se eleva al cuadrado.
# - Si los $x$ están lejos del origen (años como 2003, por ejemplo), las columnas $1,x,x^2,\dots$ son casi paralelas: $A$ está *mal condicionada* y los coeficientes en potencias de $x$ son enormes y casi se cancelan. Conviene usar $t=x-c$ con $c$ un valor central; la calculadora lo avisa y da también esa forma.
# - Las bibliotecas no forman $A^TA$, porque su número de condición es el cuadrado del de $A$ y eso empeora el redondeo: `np.linalg.lstsq` resuelve con la descomposición en valores singulares y otras usan QR. `np.polyfit` devuelve los coeficientes **de la potencia más alta a la más baja**, al revés que $\beta_0,\beta_1,\dots$

# %%
_MODELOS = ["recta", "recta por el origen", "polinomio", "exponencial", "potencia"]
_MAX_EXACTO = 12                 # con hasta 12 datos (y grado hasta 3) el cálculo es exacto
_COND_ALTA = 1e8                 # número de condición de A a partir del cual se avisa


class _Desborde(ValueError):
    """La predicción no cabe en los números decimales."""


def _modelo(modelo):
    """El nombre del modelo (acepta mayúsculas y acentos): uno de _MODELOS."""
    t = normaliza_texto(modelo) if isinstance(modelo, str) else ""
    for m in _MODELOS:
        if t == normaliza_texto(m):
            return m
    raise ValueError(f"el modelo debe ser uno de: {', '.join(_MODELOS)} (escribiste «{modelo}»)")


def _grado(grado):
    """El grado del polinomio: un entero de 1 a 6."""
    if grado is None:
        raise ValueError("indica el grado del polinomio, un entero de 1 a 6")
    g = numero(grado, "el grado")
    if not g.is_integer or not 1 <= g <= 6:
        raise ValueError(f"el grado del polinomio debe ser un entero de 1 a 6 (escribiste {_txt(g)})")
    return int(g)


def _potencia(letra, p):
    """'' (p = 0), 'x', 'x²', 'x³', ...: la misma escritura en la forma del modelo y en la ecuación."""
    return "" if p == 0 else (letra if p == 1 else f"{letra}{str(p).translate(_SUP)}")


def _forma(m, grado=1):
    """La forma del modelo con sus coeficientes: 'y = β0 + β1 x + β2 x²', 'z = β0 + β1 x + β2 y', ..."""
    if m == "recta":
        return "y = β0 + β1 x"
    if m == "recta por el origen":
        return "y = β1 x"
    if m == "polinomio":
        return "y = " + " + ".join(f"β{p}" + (f" {_potencia('x', p)}" if p else "") for p in range(grado + 1))
    if m == "exponencial":
        return "y = a·e^(k x), linealizado como ln y = ln a + k x"
    if m == "potencia":
        return "y = a·x^b, linealizado como ln y = ln a + b ln x"
    return "z = β0 + β1 x + β2 y"


def _exacto_en(c, texto):
    """True si el texto (un decimal) es exactamente el racional c, es decir, si se muestra sin redondear."""
    try:
        return isinstance(c, sp.Rational) and Fraction(texto) == Fraction(int(c.p), int(c.q))
    except (ValueError, ZeroDivisionError):
        return False


def _coef_txt(c, n):
    """Un coeficiente: '10.02' (decimal exacto), '31/14 ≈ 2.214' (fracción corta) o '≈ 9.941' (redondeado)."""
    if isinstance(c, sp.Rational):
        t = _num(c, n)
        if _exacto_en(c, t):
            return t
        if len(_txt(c)) <= 20:
            return f"{_txt(c)} ≈ {cifras(c, n)}"
    return f"≈ {cifras(c, n)}"


def _medida_txt(c, n):
    """SSE, RMSE o R²: el decimal exacto si es corto ('0.099'); si no, redondeado ('≈ 0.1407')."""
    t = _num(c, n) if isinstance(c, sp.Rational) else None
    return t if t is not None and _exacto_en(c, t) else f"≈ {cifras(c, n)}"


def _decimales_escritos(textos):
    """El mayor número de decimales con que se escribieron los datos ('10.0' -> 1, '3' -> 0), o None si alguno no
    es un decimal simple (fracción o notación científica)."""
    d = 0
    for t in textos:
        m = re.fullmatch(r"[+-]?\d*(?:\.(\d*))?", str(t).strip())
        if not m:
            return None
        d = max(d, len(m.group(1) or ""))
    return d


def _fijo(v, d):
    """v (racional exacto o decimal) redondeado (5 sube) a d decimales fijos: (14, 1) -> '14.0'; sin '-0.00'."""
    with localcontext() as ctx:
        ctx.prec = 60
        x = (Decimal(int(v.p)) / Decimal(int(v.q))) if isinstance(v, sp.Rational) else Decimal(repr(float(v)))
        t = format(x.quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP), "f")
    return t[1:] if t.startswith("-") and float(t) == 0 else t


def _decimales_exactos(v):
    """Decimales de un racional con decimal finito (3/4 -> 2), o None si no es finito (1/3) o no es racional."""
    if not isinstance(v, sp.Rational):
        return None
    q, a, b = int(v.q), 0, 0
    while q % 2 == 0:
        q, a = q // 2, a + 1
    while q % 5 == 0:
        q, b = q // 5, b + 1
    return max(a, b) if q == 1 else None


def _columna(valores, n):
    """Textos de una columna calculada (predichos o residuos), todos con los mismos decimales: si todos son decimales
    exactos que caben, tal cual (10.02, 14.00, -0.20); si no, redondeados a n cifras significativas del mayor. Con
    valores muy grandes o muy chicos, cada uno con n cifras y notación ×10^."""
    grande = max(abs(float(v)) for v in valores)
    if grande == 0:
        return ["0"] * len(valores)
    d_sig = n - 1 - math.floor(math.log10(grande))
    exactos = [_decimales_exactos(v) for v in valores]
    if all(e is not None for e in exactos) and max(exactos) <= max(d_sig, 0) and grande < 1e15:
        return [_fijo(v, max(exactos)) for v in valores]
    if grande >= 10 ** n or grande < 1e-4:
        return [cifras(v, n) for v in valores]
    return [_fijo(v, d_sig) for v in valores]


def _como_datos(v, textos, n):
    """Un valor calculado con los decimales de los datos escritos (ȳ = 14 con datos '10.0' -> '14.0'), si cabe
    exacto; si no, redondeado a n cifras. Devuelve (texto, exacto)."""
    d = _decimales_escritos(textos)
    e = _decimales_exactos(v)
    if d is not None and e is not None and e <= d:
        return _fijo(v, d), True
    t = _num(v, n)
    return t, _exacto_en(v, t)


def _redondea(c, m):
    """c redondeado (5 sube) a m cifras significativas, como número decimal."""
    with localcontext() as ctx:
        ctx.prec = 60
        if isinstance(c, sp.Rational):
            if c == 0:
                return 0.0
            d = Decimal(int(c.p)) / Decimal(int(c.q))
        else:
            v = float(c)
            if v == 0 or not math.isfinite(v):
                return v
            d = Decimal(repr(v))
        return float(d.quantize(Decimal(1).scaleb(d.adjusted() - m + 1), rounding=ROUND_HALF_UP))


def _busca_cifras(n, evalua, coefs, pred):
    """Las cifras (de n a 15) con las que una ecuación reproduce el modelo: con sus coeficientes redondeados a esas
    cifras, los predichos cambian menos que 10^-n veces el mayor de ellos. Devuelve (cifras, reproduce)."""
    pred = np.asarray(pred, dtype=float)
    escala = max(float(np.max(np.abs(pred))), 1e-300)
    for m in range(n, 16):
        with np.errstate(all="ignore"):
            v = np.asarray(evalua([_redondea(c, m) for c in coefs]), dtype=float)
        if np.all(np.isfinite(v)) and float(np.max(np.abs(v - pred))) <= 10.0 ** (-n) * escala:
            return m, True
    return n, False


def _ecuacion_lineal(izq, coef, monomios, m):
    """'y = 10.02 + 0.0199 x', 'y = 1 + x²', 'z ≈ 1.02 + 0.55 x + 1.95 y': coeficientes con m cifras; '≈' si alguno
    se redondeó. Un cero exacto no se escribe."""
    exactos, partes = True, []
    for c, mono in zip(coef, monomios):
        if isinstance(c, sp.Rational) and c == 0:
            continue
        negativo = float(c) < 0
        mag = -c if negativo else c
        t = _num(mag, m)
        exactos &= _exacto_en(mag, t)
        cuerpo = mono if (mono and t == "1") else (f"{t} {mono}" if mono else t)
        partes.append((("- " if partes else "-") if negativo else ("+ " if partes else "")) + cuerpo)
    return f"{izq} {'=' if exactos else '≈'} {' '.join(partes) if partes else '0'}"


def _ecuacion_log(m_modelo, a, k, L, m, letra="x"):
    """Exponencial y ≈ a·e^(k x) o potencia y ≈ a·x^b; si |ln a| > 700 (a no cabe en los decimales), en la forma
    y ≈ e^(ln a + k x) o y ≈ e^(ln a + b ln x)."""
    if abs(L) > 700:
        resto = f"{_num(abs(k), m)} {letra}" if m_modelo == "exponencial" else f"{_num(abs(k), m)} ln {letra}"
        return f"y ≈ e^({_num(L, m)} {'-' if k < 0 else '+'} {resto})"
    if m_modelo == "exponencial":
        return f"y ≈ {_num(a, m)}·e^({_num(k, m)} {letra})"
    return f"y ≈ {_num(a, m)}·{letra}^{_num(k, m)}"


def _a_potencias(d, c, h):
    """Coeficientes en potencias de x del polinomio d0 + d1·s + ... + dg·s^g con s = (x - c)/h. La conversión se
    hace con fracciones exactas a partir de los decimales de d, c y h; solo sirve para MOSTRAR los coeficientes:
    para evaluar el modelo se usa la forma en s, que no pierde cifras."""
    C, H = Fraction(c), Fraction(h)
    b = [Fraction(0)] * len(d)
    for p, dp in enumerate(d):
        D = Fraction(float(dp)) / H ** p
        for k in range(p + 1):
            b[k] += D * math.comb(p, k) * (-C) ** (p - k)
    return np.array([float(v) for v in b])


def _fraccion(c):
    return Fraction(int(c.p), int(c.q)) if isinstance(c, sp.Rational) else Fraction(c)


def _centro(valores, textos):
    """Un valor central 'redondo' para t = x - c: el promedio, redondeado a los decimales de los datos escritos."""
    media = sum(valores) / len(valores)
    d = _decimales_escritos(textos)
    with localcontext() as ctx:
        ctx.prec = 60
        x = Decimal(int(media.p)) / Decimal(int(media.q))
        q = x.quantize(Decimal(1).scaleb(-(d if d is not None else 4)), rounding=ROUND_HALF_UP)
    return sp.Rational(str(q))


def _medidas(y, pred, n, con_ordenada, exacto):
    """SSE, RMSE y R² (None sin ordenada al origen o si todos los y son iguales), exactos o en decimales."""
    if exacto:
        r = [a - b for a, b in zip(y, pred)]
        sse = sum(e ** 2 for e in r)
        media = sum(y) / n
        sst = sum((v - media) ** 2 for v in y)
        return sse, sp.sqrt(sse / n), (1 - sse / sst if con_ordenada and sst != 0 else None), sst
    r = np.asarray(y, float) - np.asarray(pred, float)
    sse = float(r @ r)
    yy = np.asarray(y, float)
    sst = float(((yy - yy.mean()) ** 2).sum())
    return sse, math.sqrt(sse / n), (1 - sse / sst if con_ordenada and sst > 0 else None), sst


def _revisa_finito(*arreglos):
    for a in arreglos:
        if not np.all(np.isfinite(np.asarray(a, dtype=float))):
            raise ValueError("los datos dan valores fuera del intervalo que maneja esta calculadora (desbordamiento): "
                             "usa números más pequeños o cambia de modelo")


def _aviso_cond(m, cond, n_cifras, centrada):
    """El aviso de A mal condicionada, con el consejo que corresponde al modelo; centrada dice si abajo va la forma
    centrada."""
    inicio = f"la matriz A está mal condicionada (número de condición ≈ {cifras(cond, n_cifras)}): "
    if m == "potencia":
        return (inicio + "los valores de ln x están muy juntos; cambia la unidad o la escala de x (por ejemplo, x/1000) "
                "o usa otro modelo")
    texto = inicio + ("los coeficientes son muy sensibles al redondeo" if m == "plano" else
                      "los coeficientes en potencias de x son muy sensibles al redondeo y casi se cancelan")
    if centrada:
        return texto + ("; abajo va también la forma centrada, que es la que conviene usar" if m == "plano" else
                        "; abajo va también la forma con t = x - c, con c un valor central, que es la que conviene usar")
    return texto + ("; conviene restar a x y a y valores centrales" if m == "plano" else
                    "; conviene restar a x un valor central (t = x - c)")


def ajuste(x, y, modelo, grado=None, n_cifras=4):
    """Ajuste por mínimos cuadrados de y contra x. modelo: 'recta', 'recta por el origen', 'polinomio' (con grado de
    1 a 6), 'exponencial' (y = a·e^(k x), exige y > 0) o 'potencia' (y = a·x^b, exige x > 0 y y > 0); el grado solo se
    usa con 'polinomio'. Con datos racionales, hasta 12 datos y grado hasta 3 (recta, recta por el origen o polinomio)
    el cálculo es EXACTO con las ecuaciones normales AᵀA β = Aᵀy; si no, en decimales con np.linalg.lstsq (con x
    centrado y escalado). Diccionario con A, ATA, ATy, coeficientes (β0, β1, ... o a y k / a y b), nombres, ecuacion
    (con las cifras necesarias para que reproduzca el modelo, al menos n_cifras), ecuacion_centrada (None, o la forma
    con t = x - c si A está mal condicionada), cond (número de condición de A), predichos, residuos, SSE, RMSE, R2
    (None sin ordenada al origen), avisos, exacto, modelo, x, y (y sus textos x_txt, y_txt), n y potencias. Para
    'exponencial' y 'potencia', A, ATA y ATy son los del problema en ln y, y R2 se calcula para ln y. Si las columnas
    de A son dependientes (todos los x iguales, o grado >= número de x distintos), ValueError «las columnas de A son
    dependientes: ...»."""
    m = _modelo(modelo)
    xs, x_txt = datos(x, "los datos de x", con_texto=True)
    ys, y_txt = datos(y, "los datos de y", con_texto=True)
    if len(xs) != len(ys):
        raise ValueError(f"hay {_cuenta(len(xs), 'dato', 'datos')} de x y {len(ys)} de y: debe haber un y por cada x"
                         f"{pista_coma(x) or pista_coma(y)}")
    n = len(xs)
    potencias = list(range(_grado(grado) + 1)) if m == "polinomio" else ([1] if m == "recta por el origen" else [0, 1])
    if m in ("exponencial", "potencia"):
        for lista, letra in ((ys, "y"), (xs, "x")) if m == "potencia" else ((ys, "y"),):
            malos = [j for j, v in enumerate(lista, 1) if v <= 0]
            if malos:
                forma = "y = a·e^(k x)" if m == "exponencial" else "y = a·x^b"
                raise ValueError(f"el modelo {m} {forma} necesita todos los {letra} positivos, porque se toma ln {letra}: "
                                 f"el dato {malos[0]} tiene {letra} = {(x_txt if letra == 'x' else y_txt)[malos[0] - 1]}")
    distintos = len(set(xs))
    if m == "recta por el origen":
        if all(v == 0 for v in xs):
            raise ValueError("las columnas de A son dependientes: todos los x son 0, así que la única columna de A es cero")
    elif distintos < len(potencias):
        detalle = ("todos los x son iguales, así que la columna de x es un múltiplo de la columna de unos" if distintos == 1
                   else f"con {distintos} valores distintos de x, el grado puede ser a lo más {distintos - 1}")
        raise ValueError("las columnas de A son dependientes: " + detalle)
    avisos = []
    if pista_coma(x) or pista_coma(y):
        avisos.append("hay comas entre dígitos (1,5) y números separados solo por espacios; si escribiste decimales con "
                      "coma, usa punto decimal (1.5, no 1,5): la calculadora leyó cada coma como separador")
    if len(potencias) == n:
        avisos.append({"polinomio": "el polinomio", "recta": "la recta"}.get(m, "el modelo") +
                      " pasa por todos los datos; residuo cero no significa mejor modelo")
    con_ordenada = 0 in potencias
    exacto = m in ("recta", "recta por el origen", "polinomio") and n <= _MAX_EXACTO and max(potencias) <= 3
    coef_ln, escala = None, None
    X = np.array([float(v) for v in xs])
    if exacto:
        A = sp.Matrix([[xi ** p for p in potencias] for xi in xs])
        Y = sp.Matrix(ys)
        ATA, ATy = A.T * A, A.T * Y
        coef = list(ATA.LUsolve(ATy))
        predichos = list(A * sp.Matrix(coef))
        residuos = [a - b for a, b in zip(ys, predichos)]
        sse, rmse, r2, sst = _medidas(ys, predichos, n, con_ordenada, True)
        cond = float(np.linalg.cond(_flotante(A)))
    else:
        Yv = np.array([float(v) for v in ys])
        Xl = np.log(X) if m == "potencia" else X
        Yl = np.log(Yv) if m in ("exponencial", "potencia") else Yv
        with np.errstate(all="ignore"):
            A = np.column_stack([Xl ** p for p in potencias])
            _revisa_finito(A)
            ATA, ATy = A.T @ A, A.T @ Yl
            cond = float(np.linalg.cond(A))
            if con_ordenada:                      # x centrado y escalado: mismas columnas generadas, bien condicionadas
                c0 = float(np.mean(Xl))
                h = float(np.max(np.abs(Xl - c0))) or 1.0
                S = np.column_stack([((Xl - c0) / h) ** p for p in potencias])
                d = np.linalg.lstsq(S, Yl, rcond=None)[0]
                escala = (c0, h, [float(v) for v in d])
                cl = _a_potencias(d, c0, h)
                ajuste_l = S @ d
            else:
                cl = np.linalg.lstsq(A, Yl, rcond=None)[0]
                ajuste_l = A @ cl
            if m in ("exponencial", "potencia"):
                coef_ln = [float(v) for v in cl]
                coef = [float(np.exp(cl[0])), float(cl[1])]
                predichos = np.exp(ajuste_l)
                avisos.append("el ajuste se hace con ln y (una recta en ln y): minimiza el error en ln y, no en y; "
                              "R² se calcula para ln y")
                if abs(coef_ln[0]) > 700:
                    consejo = ("mide x desde otro origen (abajo va la forma con t = x - c)" if m == "exponencial" else
                               "cambia la unidad o la escala de x (por ejemplo, x/1000) o usa otro modelo")
                    avisos.append(f"a = e^(ln a), con ln a ≈ {cifras(coef_ln[0], n_cifras)}, no cabe en los números "
                                  f"decimales (o casi): {consejo}")
                _revisa_finito(coef[1:], predichos, ATA, ATy)
            else:
                coef = [float(v) for v in cl]
                predichos = ajuste_l
                _revisa_finito(coef, predichos, ATA, ATy)
        predichos = [float(v) for v in predichos]
        residuos = [float(a) - b for a, b in zip(ys, predichos)]
        sse, rmse, r2, sst = _medidas(ys, predichos, n, con_ordenada, False)
        if m in ("exponencial", "potencia"):
            _, _, r2, _ = _medidas(list(Yl), list(ajuste_l), n, True, False)
    if con_ordenada and (sst == 0) and m not in ("exponencial", "potencia"):
        avisos.append("todos los y son iguales: R² no está definido")
    # ecuación en x, con las cifras que hagan falta para reproducir el modelo
    pf = [float(v) for v in predichos]
    if m in ("exponencial", "potencia"):
        lx = np.log(X) if m == "potencia" else X
        if abs(coef_ln[0]) > 700:
            cifras_ec, reproduce = _busca_cifras(n_cifras, lambda c: np.exp(c[0] + c[1] * lx), coef_ln, pf)
        else:
            cifras_ec, reproduce = _busca_cifras(n_cifras, lambda c: c[0] * np.exp(c[1] * lx), coef, pf)
        ecuacion = _ecuacion_log(m, coef[0], coef[1], coef_ln[0], cifras_ec)
    else:
        cifras_ec, reproduce = _busca_cifras(n_cifras, lambda c: sum(ci * X ** p for ci, p in zip(c, potencias)), coef, pf)
        ecuacion = _ecuacion_lineal("y", coef, [_potencia("x", p) for p in potencias], cifras_ec)
    # forma centrada t = x - c, si se lee mejor que la de x
    ecuacion_centrada = None
    if m != "potencia" and con_ordenada and (cond > _COND_ALTA or cifras_ec > n_cifras or not reproduce
                                            or (coef_ln is not None and abs(coef_ln[0]) > 700)):
        c = _centro(xs, x_txt)
        T = X - float(c)
        if m == "exponencial":
            c0, h, d = escala
            Lc, kc = d[0] + d[1] * (float(c) - c0) / h, d[1] / h
            cc, rep = _busca_cifras(n_cifras, lambda q: q[0] * np.exp(q[1] * T), [math.exp(Lc), kc], pf)
            if rep and (cc < cifras_ec or not reproduce or abs(coef_ln[0]) > 700):
                ecuacion_centrada = _ecuacion_log(m, math.exp(Lc), kc, Lc, cc, "t") + f", con t = x - {_num(c, 15)}"
        else:
            if exacto:
                gam = [sum(coef[k] * sp.binomial(k, j) * c ** (k - j) for k in range(j, len(potencias)))
                       for j in range(len(potencias))]
            else:
                c0, h, d = escala
                gam = list(_a_potencias(d, Fraction(c0) - _fraccion(c), h))
            cc, rep = _busca_cifras(n_cifras, lambda q: sum(qi * T ** p for p, qi in enumerate(q)), gam, pf)
            if rep and (cc < cifras_ec or not reproduce):
                ecuacion_centrada = (_ecuacion_lineal("y", gam, [_potencia("t", p) for p in potencias], cc) +
                                     f", con t = x - {_num(c, 15)}")
    if cond > _COND_ALTA:
        avisos.append(_aviso_cond(m, cond, n_cifras, ecuacion_centrada is not None))
    if cifras_ec > n_cifras:
        avisos.append(f"la ecuación se escribe con {cifras_ec} cifras, y no con {n_cifras}, para que reproduzca el "
                      "modelo: sus términos casi se cancelan")
    if not reproduce:
        avisos.append("con los coeficientes redondeados, la ecuación en x no reproduce el modelo (sus términos casi se "
                      "cancelan)" + (": usa la forma centrada" if ecuacion_centrada else ""))
    nombres = {"exponencial": ["a", "k"], "potencia": ["a", "b"]}.get(m, [f"β{p}" for p in potencias])
    return {"modelo": m, "grado": max(potencias), "x": xs, "y": ys, "x_txt": x_txt, "y_txt": y_txt, "n": n, "A": A,
            "ATA": ATA, "ATy": ATy, "coeficientes": coef, "coeficientes_ln": coef_ln, "nombres": nombres,
            "ecuacion": ecuacion, "cifras_ecuacion": cifras_ec, "ecuacion_reproduce": reproduce,
            "ecuacion_centrada": ecuacion_centrada, "cond": cond, "predichos": predichos, "residuos": residuos,
            "SSE": sse, "RMSE": rmse, "R2": r2, "avisos": avisos, "exacto": exacto, "potencias": potencias,
            "escala": escala}


def ajuste_plano(x, y, z, n_cifras=4):
    """Ajuste del plano z = β0 + β1 x + β2 y por mínimos cuadrados (exacto con hasta 12 datos; con más, en decimales
    con x y y centrados y escalados). El mismo diccionario que ajuste(), con modelo 'plano', y_datos (los y) y z.
    Si los puntos (x, y) están sobre una recta, ValueError."""
    xs, x_txt = datos(x, "los datos de x", con_texto=True)
    ys, y_txt = datos(y, "los datos de y", con_texto=True)
    zs, z_txt = datos(z, "los datos de z", con_texto=True)
    if not len(xs) == len(ys) == len(zs):
        raise ValueError(f"hay {_cuenta(len(xs), 'dato', 'datos')} de x, {len(ys)} de y y {len(zs)} de z: deben ser "
                         "los mismos")
    n = len(xs)
    G = sp.Matrix([[n, sum(xs), sum(ys)], [sum(xs), sum(a * a for a in xs), sum(a * b for a, b in zip(xs, ys))],
                   [sum(ys), sum(a * b for a, b in zip(xs, ys)), sum(b * b for b in ys)]])
    if G.det() == 0:
        raise ValueError("las columnas de A son dependientes: los puntos (x, y) están todos sobre una misma recta (o son "
                         "el mismo punto), así que el plano no queda determinado")
    avisos = ["el plano pasa por los tres datos; residuo cero no significa mejor modelo"] if n == 3 else []
    exacto = n <= _MAX_EXACTO
    X, Yd = np.array([float(v) for v in xs]), np.array([float(v) for v in ys])
    escala = None
    if exacto:
        A = sp.Matrix([[1, a, b] for a, b in zip(xs, ys)])
        Z = sp.Matrix(zs)
        ATA, ATy = A.T * A, A.T * Z
        coef = list(ATA.LUsolve(ATy))
        predichos = list(A * sp.Matrix(coef))
        residuos = [a - b for a, b in zip(zs, predichos)]
        cond = float(np.linalg.cond(_flotante(A)))
    else:
        A = np.column_stack([np.ones(n), X, Yd])
        Zf = np.array([float(v) for v in zs])
        ATA, ATy = A.T @ A, A.T @ Zf
        cond = float(np.linalg.cond(A))
        cx, cy = float(X.mean()), float(Yd.mean())
        hx, hy = float(np.max(np.abs(X - cx))) or 1.0, float(np.max(np.abs(Yd - cy))) or 1.0
        S = np.column_stack([np.ones(n), (X - cx) / hx, (Yd - cy) / hy])
        d = np.linalg.lstsq(S, Zf, rcond=None)[0]
        b1, b2 = Fraction(float(d[1])) / Fraction(hx), Fraction(float(d[2])) / Fraction(hy)
        coef = [float(Fraction(float(d[0])) - b1 * Fraction(cx) - b2 * Fraction(cy)), float(b1), float(b2)]
        predichos = [float(v) for v in S @ d]
        escala = (cx, hx, cy, hy, [float(v) for v in d])
        _revisa_finito(coef, predichos)
        residuos = [float(a) - b for a, b in zip(zs, predichos)]
    sse, rmse, r2, sst = _medidas(zs, predichos, n, True, exacto)
    if sst == 0:
        avisos.append("todos los z son iguales: R² no está definido")
    Xc0, Yc0 = X - X.mean(), Yd - Yd.mean()
    colineal = float(np.linalg.cond(np.column_stack([Xc0 / np.linalg.norm(Xc0), Yc0 / np.linalg.norm(Yc0)])))
    if colineal > 1e3:
        avisos.append("los puntos (x, y) están casi sobre una recta: el plano queda mal determinado y sus coeficientes "
                      "son muy sensibles a los datos (un cambio pequeño en un dato los cambia mucho)")
    pf = [float(v) for v in predichos]
    cifras_ec, reproduce = _busca_cifras(n_cifras, lambda c: c[0] + c[1] * X + c[2] * Yd, coef, pf)
    ecuacion_centrada = None
    if cond > _COND_ALTA or cifras_ec > n_cifras or not reproduce:
        cx_, cy_ = _centro(xs, x_txt), _centro(ys, y_txt)
        g0 = (coef[0] + coef[1] * cx_ + coef[2] * cy_) if exacto else float(
            Fraction(coef[0]) + Fraction(coef[1]) * _fraccion(cx_) + Fraction(coef[2]) * _fraccion(cy_))
        gam = [g0, coef[1], coef[2]]
        Xc, Yc = X - float(cx_), Yd - float(cy_)
        cc, rep = _busca_cifras(n_cifras, lambda q: q[0] + q[1] * Xc + q[2] * Yc, gam, pf)
        if rep and (cc < cifras_ec or not reproduce):
            ecuacion_centrada = _ecuacion_lineal("z", gam, ["", f"(x - {_num(cx_, 15)})", f"(y - {_num(cy_, 15)})"], cc)
    if cond > _COND_ALTA:
        avisos.append(_aviso_cond("plano", cond, n_cifras, ecuacion_centrada is not None))
    if cifras_ec > n_cifras:
        avisos.append(f"la ecuación se escribe con {cifras_ec} cifras, y no con {n_cifras}, para que reproduzca el "
                      "modelo: sus términos casi se cancelan")
    if not reproduce:
        avisos.append("con los coeficientes redondeados, la ecuación no reproduce el modelo"
                      + (": usa la forma centrada" if ecuacion_centrada else ""))
    return {"modelo": "plano", "grado": 1, "x": xs, "y_datos": ys, "y": zs, "z": zs, "x_txt": x_txt, "yd_txt": y_txt,
            "y_txt": z_txt, "n": n, "A": A, "ATA": ATA, "ATy": ATy, "coeficientes": coef, "coeficientes_ln": None,
            "nombres": ["β0", "β1", "β2"], "ecuacion": _ecuacion_lineal("z", coef, ["", "x", "y"], cifras_ec),
            "cifras_ecuacion": cifras_ec, "ecuacion_reproduce": reproduce, "ecuacion_centrada": ecuacion_centrada,
            "cond": cond, "predichos": predichos, "residuos": residuos, "SSE": sse, "RMSE": rmse, "R2": r2,
            "avisos": avisos, "exacto": exacto, "potencias": None, "escala": escala}


def _evalua_flotante(res, t):
    """El modelo evaluado en un arreglo t (para dibujar y predecir en decimales); lo que se desborda queda como nan."""
    m, c = res["modelo"], [float(v) for v in res["coeficientes"]]
    t = np.asarray(t, dtype=float)
    with np.errstate(all="ignore"):
        tl = np.log(t) if m == "potencia" else t
        if res.get("escala"):                       # la forma en s = (x - c)/h, sin pérdida de cifras
            c0, h, d = res["escala"]
            s = (tl - c0) / h
            v = sum(dp * s ** p for p, dp in enumerate(d))
        elif m in ("exponencial", "potencia"):
            v = res["coeficientes_ln"][0] + c[1] * tl
        else:
            v = sum(ci * tl ** p for ci, p in zip(c, res["potencias"]))
        if m in ("exponencial", "potencia"):
            v = np.exp(v)
    v = np.array(v, dtype=float)
    v[~np.isfinite(v)] = np.nan
    return v


def _evalua(res, x0):
    """El modelo evaluado en x0: exacto si el ajuste fue exacto; si no, en decimales (_Desborde si no cabe)."""
    if res["exacto"]:
        return sum(ci * x0 ** p for ci, p in zip(res["coeficientes"], res["potencias"]))
    v = float(_evalua_flotante(res, np.array([float(x0)]))[0])
    if not math.isfinite(v):
        raise _Desborde(f"la predicción en x0 = {_num(x0, 15)} se desborda: el valor no cabe en los números decimales")
    return v


def predice(resultado, x0):
    """Predicción del modelo en x0: diccionario con x0, y (exacto si el ajuste fue exacto) y aviso (None, o
    'extrapolación: x0 = ... está fuera del intervalo de los datos [min, max]...')."""
    if resultado["modelo"] == "plano":
        raise ValueError("predice es para los modelos de una variable; para el plano evalúa z = β0 + β1 x + β2 y")
    x0 = numero(x0, "x0")
    if resultado["modelo"] == "potencia" and x0 <= 0:
        raise ValueError("el modelo potencia y = a·x^b solo se evalúa en x positivos")
    lo, hi = min(resultado["x"]), max(resultado["x"])
    aviso = None
    if not lo <= x0 <= hi:
        aviso = (f"extrapolación: x0 = {_num(x0, 15)} está fuera del intervalo de los datos [{_num(lo, 15)}, "
                 f"{_num(hi, 15)}]; nada garantiza que el modelo valga ahí")
    return {"x0": x0, "y": _evalua(resultado, x0), "aviso": aviso}


def compara_bibliotecas(x, y, grado):
    """Coeficientes de un polinomio de grado dado por tres caminos: ecuaciones normales (np.linalg.solve con AᵀA),
    np.linalg.lstsq (descomposición en valores singulares, sin AᵀA) y np.polyfit (que los da AL REVÉS, de la
    potencia más alta a la más baja). Diccionario con propios (los de ajuste()), normales, lstsq, polyfit,
    cond_A y cond_ATA (los números de condición de A y de AᵀA: el segundo es el cuadrado del primero), coinciden y
    texto."""
    res = ajuste(x, y, "polinomio", grado)
    g = res["grado"]
    X, Y = np.array([float(v) for v in res["x"]]), np.array([float(v) for v in res["y"]])
    A = np.column_stack([X ** p for p in range(g + 1)])
    cond_A, cond_ATA = float(np.linalg.cond(A)), float(np.linalg.cond(A.T @ A))
    try:
        normales = np.linalg.solve(A.T @ A, A.T @ Y)
    except np.linalg.LinAlgError:
        raise ValueError("AᵀA resultó singular en decimales: A está demasiado mal condicionada para las ecuaciones "
                         "normales") from None
    lstsq = np.linalg.lstsq(A, Y, rcond=None)[0]
    polyfit = np.polyfit(X, Y, g)
    propios = np.array([float(c) for c in res["coeficientes"]])
    escala = max(1.0, float(np.abs(propios).max()))
    coinciden = all(np.allclose(v, propios, rtol=0, atol=1e-9 * escala) for v in (normales, lstsq, polyfit[::-1]))
    orden = ", ".join(f"β{p}" for p in range(g, -1, -1))
    al_reves = ", ".join(f"β{p}" for p in range(g + 1))
    texto = (f"np.polyfit devuelve ({orden}): de la potencia más alta a la más baja. Leídos al revés, "
             f"({al_reves}), son los mismos que dan las ecuaciones normales y np.linalg.lstsq."
             + ("" if coinciden else " (Aquí difieren más de lo esperado: A está mal condicionada.)")
             + f" Números de condición: cond(A) ≈ {cifras(cond_A, 4)} y cond(AᵀA) ≈ {cifras(cond_ATA, 4)}, que es "
               "cond(A)²: por eso las bibliotecas no forman AᵀA.")
    return {"propios": propios, "normales": normales, "lstsq": lstsq, "polyfit": polyfit, "cond_A": cond_A,
            "cond_ATA": cond_ATA, "coinciden": coinciden, "texto": texto}


# --- informes ---------------------------------------------------------------------------------------
def _num_exacto(e, n):
    """Una entrada exacta de AᵀA o de Aᵀy: tal cual si su decimal es finito y corto (12015, 159.9); si no, con n
    cifras y notación ×10^ (nunca un entero redondeado que parezca exacto)."""
    d = _decimales_exactos(e)
    if d is not None and len(_fijo(e, d).replace("-", "").replace(".", "").lstrip("0")) <= 20:
        return _fijo(e, d)
    return cifras(e, n)


def _muestra_num(nombre, M, n, tex, referencia=False):
    """Muestra AᵀA o Aᵀy: exactas (= si ninguna entrada se redondeó, ≈ si alguna sí) o, si son de numpy,
    redondeadas a n cifras."""
    if isinstance(M, np.ndarray):
        muestra_decimales(nombre, M, n, tex)
        return
    textos = [[_num_exacto(M[i, j], n) for j in range(M.cols)] for i in range(M.rows)]
    igual = all(_exacto_en(M[i, j], textos[i][j]) for i in range(M.rows) for j in range(M.cols))
    signo, signo_tex = ("=", "=") if igual else ("≈", r"\approx")
    filas_tex = r" \\ ".join(" & ".join(_tex_cifra(t) for t in fila) for fila in textos)
    anchos = [max(len(textos[i][j]) for i in range(M.rows)) for j in range(M.cols)]
    plano = "\n".join("[ " + "  ".join(t.rjust(a) for t, a in zip(fila, anchos)) + " ]" for fila in textos)
    _emite(_Formula(f"{tex} {signo_tex} " + r"\left[\begin{matrix}" + filas_tex + r"\end{matrix}\right]" +
                    f"\\qquad ({M.rows}\\times {M.cols})", f"{nombre}  ({dims(M)}) {signo}\n{plano}"))


def _tabla(res, n_cifras, max_filas=20):
    """Tabla de datos (como se escribieron), predichos y residuos (con los mismos decimales en cada columna)."""
    plano = res["modelo"] == "plano"
    letra, sombrero = ("z", "ẑ") if plano else ("y", "ŷ")
    k = min(res["n"], max_filas)
    columnas = [[str(i + 1) for i in range(k)], res["x_txt"][:k]] + ([res["yd_txt"][:k]] if plano else [])
    columnas += [res["y_txt"][:k], _columna(res["predichos"][:k], n_cifras), _columna(res["residuos"][:k], n_cifras)]
    enc = ["dato", "x"] + (["y"] if plano else []) + [letra, f"{sombrero} (predicho)", f"r = {letra} - {sombrero}"]
    anchos = [max(len(t) for t in col + [e]) for col, e in zip(columnas, enc)]
    print("  ".join(e.rjust(a) for e, a in zip(enc, anchos)))
    for i in range(k):
        print("  ".join(col[i].rjust(a) for col, a in zip(columnas, anchos)))
    if res["n"] > max_filas:
        print(f"  ... ({res['n']} datos; la tabla muestra los primeros {max_filas})")


def _informe_ajuste(res, n_cifras, prediccion=None):
    m, n = res["modelo"], res["n"]
    letra = "z" if m == "plano" else "y"
    nombre = f"polinomio de grado {res['grado']}" if m == "polinomio" else m
    print(f"Modelo: {nombre}, {_forma(m, res['grado'])}; {_cuenta(n, 'dato', 'datos')}.")
    if res["exacto"]:
        print(f"Cálculo exacto con fracciones: ecuaciones normales AᵀA β = Aᵀ{letra}.")
    else:
        razon = ("el modelo usa logaritmos" if m in ("exponencial", "potencia") else
                 f"hay más de {_MAX_EXACTO} datos" if n > _MAX_EXACTO else "el grado es mayor que 3")
        print(f"Cálculo en decimales con np.linalg.lstsq (sin formar AᵀA), porque {razon}.")
    cols = {"exponencial": ["1", "x"], "potencia": ["1", "ln x"], "plano": ["1", "x", "y"]}.get(
        m, [_potencia("x", p) or "1" for p in (res["potencias"] or [])])
    print(f"A tiene {_cuenta(n, 'renglón', 'renglones')} (uno por dato) y {_cuenta(len(cols), 'columna', 'columnas')}: "
          f"{_lista(cols)}" + ("; el lado derecho es ln y." if m in ("exponencial", "potencia") else "."))
    ATA = res["ATA"]
    k = ATA.rows if isinstance(ATA, sp.MatrixBase) else ATA.shape[0]
    derecha, tex_der = (("Aᵀ(ln y)", r"A^T(\ln\mathbf y)") if m in ("exponencial", "potencia") else
                        (f"Aᵀ{letra}", rf"A^T\mathbf {letra}"))
    if k <= 4:
        if not res["exacto"]:
            print("AᵀA y " + derecha + ", solo como referencia (lstsq no los usa):")
        _muestra_num("AᵀA", ATA, n_cifras, r"A^TA")
        _muestra_num(derecha, res["ATy"], n_cifras, tex_der)
    else:
        print(f"AᵀA es de {k}×{k}: no se muestra aquí (está en el resultado de ajuste, clave 'ATA').")
    if res["coeficientes_ln"] is not None:
        b = "k" if m == "exponencial" else "b"
        print(f"Recta en ln y: ln a ≈ {cifras(res['coeficientes_ln'][0], n_cifras)}, {b} ≈ "
              f"{cifras(res['coeficientes_ln'][1], n_cifras)};  a = e^(ln a)")
    partes = []
    for nm, c in zip(res["nombres"], res["coeficientes"]):
        if nm == "a" and abs(res["coeficientes_ln"][0]) > 700:
            partes.append(f"a = e^(ln a), con ln a ≈ {cifras(res['coeficientes_ln'][0], n_cifras)} (no cabe en los "
                          "números decimales)")
        else:
            partes.append(_eq(nm, _coef_txt(c, n_cifras)))
    print("Coeficientes: " + ",  ".join(partes))
    print("Ecuación: " + res["ecuacion"])
    if res["ecuacion_centrada"]:
        print("Forma centrada (se lee mejor): " + res["ecuacion_centrada"])
    _tabla(res, n_cifras)
    if m in ("recta", "polinomio", "plano"):
        suma = sum(res["residuos"])
        cero = suma == 0 if res["exacto"] else abs(suma) <= 1e-9 * max(1.0, max(abs(float(v)) for v in res["y"]))
        print(_eq("Suma de residuos", _num(suma, n_cifras) if res["exacto"] else f"≈ {cifras(suma, n_cifras)}") +
              ((": es cero por la columna de unos." if res["exacto"] else
                ": es cero salvo el redondeo, por la columna de unos.") if cero else ""))
        if m == "recta":
            tx, ex = _como_datos(sum(res["x"]) / n, res["x_txt"], n_cifras)
            ty, ey = _como_datos(sum(res["y"]) / n, res["y_txt"], n_cifras)
            print(f"La recta pasa por (x̄, ȳ) {'=' if ex and ey else '≈'} ({tx}, {ty}).")
    print(f"{_eq('SSE', _medida_txt(res['SSE'], n_cifras))},  {_eq('RMSE', _medida_txt(res['RMSE'], n_cifras))}" +
          (f",  {_eq('R²', _medida_txt(res['R2'], n_cifras))}" + (" (de ln y)" if m in ("exponencial", "potencia") else "")
           if res["R2"] is not None else ",  R²: no se da (el modelo no tiene ordenada al origen o todos los "
                                         f"{letra} son iguales)"))
    if res["exacto"]:
        print("SSE = 0: los datos están exactamente sobre el modelo (el sistema sí tenía solución)." if res["SSE"] == 0
              else "SSE > 0: ningún modelo de esta forma pasa por todos los datos; el ajuste es la solución de "
                   "mínimos cuadrados.")
    for aviso in res["avisos"]:
        print("Aviso:", aviso + ".")
    if prediccion is not None:
        print(f"Predicción en x = {_num(prediccion['x0'], 15)}: {_eq('ŷ', _coef_txt(prediccion['y'], n_cifras))}")
        if prediccion["aviso"]:
            print("Aviso:", prediccion["aviso"] + ".")


# --- dibujos ----------------------------------------------------------------------------------------
def _residuos_dibujo(res, n_cifras=4):
    """Los residuos para dibujar y el título del panel: si el modelo pasa por todos los datos en decimales, los
    residuos son solo redondeo y se dibujan como 0."""
    R = np.array([float(v) for v in res["residuos"]])
    s = _medida_txt(res["SSE"], n_cifras)
    titulo = f"Residuos: SSE {s}" if s.startswith("≈") else f"Residuos: SSE = {s}"
    escala = max(1.0, max(abs(float(v)) for v in res["y"]))
    if not res["exacto"] and np.all(np.abs(R) <= 1e-9 * escala):
        return np.zeros_like(R), "Residuos: 0 (el modelo pasa por todos los datos; lo demás es redondeo)"
    return R, titulo


def dibuja_ajuste(resultado, etiqueta_x="x", etiqueta_y="y", unidad_y="", prediccion=None, n_cifras=4):
    """Dos paneles: los datos con el modelo (curva suave, 400 puntos) y los residuos (tallos con la línea en cero).
    etiqueta_x y etiqueta_y rotulan los ejes (magnitud y unidad); unidad_y es la unidad de los residuos."""
    if resultado["modelo"] == "plano":
        return _dibuja_plano(resultado, n_cifras)
    X = np.array([float(v) for v in resultado["x"]])
    Y = np.array([float(v) for v in resultado["y"]])
    R, titulo_r = _residuos_dibujo(resultado, n_cifras)
    lo, hi = float(X.min()), float(X.max())
    if prediccion is not None:
        lo, hi = min(lo, float(prediccion["x0"])), max(hi, float(prediccion["x0"]))
    if resultado["modelo"] == "recta por el origen":
        lo, hi = min(lo, 0.0), max(hi, 0.0)
    ancho = (hi - lo) or max(1.0, abs(hi))
    a, b = lo - 0.05 * ancho, hi + 0.05 * ancho
    if resultado["modelo"] == "potencia" and a <= 0:
        a = 0.5 * lo if lo > 0 else float(X.min()) * 0.5
    t = np.linspace(a, b, 400)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.8))
    ax1.plot(t, _evalua_flotante(resultado, t), "-", color=_COLORES[0], lw=2.3,
             label=textwrap.fill(f"modelo: {resultado['ecuacion']}", 48))
    ax1.plot(X, Y, "o", color="black", ms=6, zorder=3, label="datos")
    if prediccion is not None:
        ax1.plot(float(prediccion["x0"]), float(prediccion["y"]), "X", color=_COLORES[1], ms=11, zorder=4,
                 label=f"predicción en x = {_num(prediccion['x0'], 15)}" + (" (extrapolación)" if prediccion["aviso"] else ""))
    ax1.set_xlabel(etiqueta_x); ax1.set_ylabel(etiqueta_y); ax1.grid(alpha=0.4)
    ax1.set_title(f"Datos y modelo ajustado ({_cuenta(resultado['n'], 'dato', 'datos')})", fontsize=10)
    ax1.legend(fontsize=8, loc="best")
    ax2.vlines(X, 0, R, colors="0.35", linestyles="-", lw=1.6)
    ax2.plot(X, R, "s", color=_COLORES[2], ms=6, zorder=3, label="residuo r = y - ŷ (tallo y cuadro)")
    ax2.axhline(0, color="black", lw=1.0, ls="--", label="r = 0")
    ax2.set_xlabel(etiqueta_x)
    ax2.set_ylabel("residuo y - ŷ" + (f" ({unidad_y})" if unidad_y else ""))
    ax2.grid(alpha=0.4)
    ax2.set_title(textwrap.fill(titulo_r, 55), fontsize=10)
    ax2.legend(fontsize=8, loc="best")
    plt.tight_layout()
    plt.show()


def _dibuja_plano(res, n_cifras):
    """Plano ajustado con los datos y sus residuos (segmentos verticales), y los residuos por número de dato."""
    X = np.array([float(v) for v in res["x"]]); Y = np.array([float(v) for v in res["y_datos"]])
    Z = np.array([float(v) for v in res["z"]]); Zp = np.array([float(v) for v in res["predichos"]])
    fig = plt.figure(figsize=(11.5, 5.0))
    ax = fig.add_subplot(1, 2, 1, projection="3d")
    gx = np.linspace(X.min(), X.max(), 12) if X.max() > X.min() else np.array([X.min() - 1, X.min() + 1])
    gy = np.linspace(Y.min(), Y.max(), 12) if Y.max() > Y.min() else np.array([Y.min() - 1, Y.min() + 1])
    GX, GY = np.meshgrid(gx, gy)
    b0, b1, b2 = (float(v) for v in res["coeficientes"])
    with np.errstate(all="ignore"):
        GZ = b0 + b1 * GX + b2 * GY
    ax.plot_surface(GX, GY, GZ, color="0.7", alpha=0.35, edgecolor="0.45", linewidth=0.4)
    ax.scatter(X, Y, Z, color="black", s=30, depthshade=False, label="datos")
    for xi, yi, zi, zp in zip(X, Y, Z, Zp):
        ax.plot([xi, xi], [yi, yi], [zi, zp], ":", color=_COLORES[1], lw=1.8)
    ax.plot([], [], ":", color=_COLORES[1], lw=1.8, label="residuo (vertical, punteado)")
    ax.plot([], [], "s", color="0.7", label=textwrap.fill(f"plano: {res['ecuacion']}", 48))
    ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
    ax.set_title("Datos y plano ajustado", fontsize=10)
    ax.legend(fontsize=8, loc="upper left")
    ax2 = fig.add_subplot(122)
    k = np.arange(1, res["n"] + 1)
    R, titulo_r = _residuos_dibujo(res, n_cifras)
    ax2.vlines(k, 0, R, colors="0.35", lw=1.6)
    ax2.plot(k, R, "s", color=_COLORES[2], ms=6, label="residuo r = z - ẑ")
    ax2.axhline(0, color="black", lw=1.0, ls="--", label="r = 0")
    ax2.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(integer=True))     # números de dato enteros
    ax2.set_xlabel("número de dato"); ax2.set_ylabel("residuo z - ẑ"); ax2.grid(alpha=0.4)
    ax2.set_title(textwrap.fill(titulo_r, 55), fontsize=10)
    ax2.legend(fontsize=8)
    plt.tight_layout()
    plt.show()

# %% [markdown]
# **Recorrido paso a paso.** La semilla del libro: un resorte cuelga con masas de 0, 100, 200, 300 y 400 g y mide 10.0, 12.1, 13.8, 16.2 y 17.9 cm. Los cinco puntos no están sobre una recta, así que el sistema $A\boldsymbol\beta=\mathbf y$ no tiene solución. Mínimos cuadrados da $y=10.02+0.0199\,x$ (con $x$ en gramos y $y$ en centímetros), con residuos $(-0.02,\,0.09,\,-0.20,\,0.21,\,-0.08)$ y $\mathrm{SSE}=0.099$. La recta pasa por $(\bar x,\bar y)=(200,\,14.0)$ y por ninguno de los datos; la recta que une el primer y el último dato tiene $\mathrm{SSE}=0.11375$, más grande. La pendiente da la constante del resorte: 100 g (peso $0.1\cdot 9.81=0.981$ N) alargan $1.99$ cm $=0.0199$ m, así que $k=\dfrac{0.981}{0.0199}\approx 49.3$ N/m.

# %%
_xs, _ys = "0 100 200 300 400", "10.0 12.1 13.8 16.2 17.9"
_semilla = ajuste(_xs, _ys, "recta")
_informe_ajuste(_semilla, 4)
_X, _Y = datos(_xs), datos(_ys)
_xb, _yb = sum(_X) / 5, sum(_Y) / 5
_Sxy = sum((a - _xb) * (b - _yb) for a, b in zip(_X, _Y))
_Sxx = sum((a - _xb) ** 2 for a in _X)
print(f"\nEn forma cerrada: x̄ = {_num(_xb)}, ȳ = {_fijo(_yb, 1)}, Sxy = {_num(_Sxy)}, Sxx = {_num(_Sxx)};  "
      f"β1 = Sxy/Sxx = {_num(_Sxy / _Sxx)},  β0 = ȳ - β1·x̄ = {_num(_yb - _Sxy / _Sxx * _xb)}")
print("¿La recta pasa por algún dato?", "sí" if any(r == 0 for r in _semilla["residuos"]) else "no, por ninguno")
_m = (_Y[-1] - _Y[0]) / (_X[-1] - _X[0])
_sse_extremos = sum((b - (_Y[0] + _m * (a - _X[0]))) ** 2 for a, b in zip(_X, _Y))
print(f"Recta que une el primer y el último dato: y = {_fijo(_Y[0], 1)} + {_num(_m)} x, con SSE = {_num(_sse_extremos, 6)} "
      f"> {_num(_semilla['SSE'])}")
_b1 = _semilla["coeficientes"][1]
_k = sp.Rational(981, 1000) / _b1          # 100 g alargan 100·β1 cm = β1 m (numéricamente)
print(f"Constante del resorte: 100 g (0.1 kg, peso 0.1·9.81 = 0.981 N) alargan 100·{_num(_b1)} = {_num(100 * _b1)} cm "
      f"= {_num(_b1)} m, así que k = 0.981/{_num(_b1)} ≈ {cifras(_k, 3)} N/m")
dibuja_ajuste(_semilla, "masa (g)", "longitud del resorte (cm)", "cm")

# %% [markdown]
# **Las ecuaciones normales a la vista.** Con los datos $(0,1)$, $(1,3)$, $(2,4)$ y $(3,4)$ y el modelo $y=\beta_0+\beta_1x$:
#
# $$A=\begin{bmatrix}1&0\\1&1\\1&2\\1&3\end{bmatrix},\qquad A^TA=\begin{bmatrix}4&6\\6&14\end{bmatrix},\qquad A^T\mathbf y=\begin{bmatrix}12\\23\end{bmatrix},$$
#
# y resolver $A^TA\,\boldsymbol\beta=A^T\mathbf y$ da $(\beta_0,\beta_1)=(3/2,\,1)$: $y=\tfrac32+x$. Al final, la descarga de un capacitor: $t=0,1,2,3$ s y $v=10.0,\ 6.1,\ 3.6,\ 2.3$ V. El modelo $v=V_0\,e^{kt}$ se linealiza como $\ln v=\ln V_0+kt$; sale $V_0\approx 9.94$ V y la constante de tiempo $\tau=-1/k\approx 2.03$ s. También se comparan las bibliotecas con la semilla, con los números de condición de $A$ y de $A^TA$.

# %%
_r2 = ajuste("0 1 2 3", "1 3 4 4", "recta")
_informe_ajuste(_r2, 4)
print("\n" + "=" * 90 + "\n\nDescarga de un capacitor:")
_desc = ajuste("0 1 2 3", "10.0 6.1 3.6 2.3", "exponencial")
_informe_ajuste(_desc, 4)
_V0, _kd = _desc["coeficientes"]
print(f"V0 = a ≈ {cifras(_V0, 3)} V,  k ≈ {cifras(_kd, 3)} 1/s,  τ = -1/k ≈ {cifras(-1 / _kd, 3)} s")
dibuja_ajuste(_desc, "tiempo t (s)", "voltaje v (V)", "V")
print("\n" + "=" * 90 + "\n\nLa semilla con tres caminos de numpy:")
_cmp = compara_bibliotecas(_xs, _ys, 1)
muestra_decimales("ecuaciones normales (β0, β1)", _cmp["normales"], 6, tex=r"\text{ecuaciones normales }(\beta_0,\ \beta_1)")
muestra_decimales("np.linalg.lstsq (β0, β1)", _cmp["lstsq"], 6, tex=r"\text{np.linalg.lstsq }(\beta_0,\ \beta_1)")
muestra_decimales("np.polyfit (β1, β0)", _cmp["polyfit"], 6, tex=r"\text{np.polyfit }(\beta_1,\ \beta_0)")
print(_cmp["texto"])

# %% [markdown]
# **Calculadora de ajuste.** Escribe los $x$ y los $y$ como listas de números separados por espacios o comas (el mismo número de datos), elige el modelo y, para el polinomio, el grado (el control de grado solo se usa con «polinomio»). Si escribes un valor en «predecir en x =», la calculadora lo evalúa y avisa si es una extrapolación. Con hasta 12 datos y grado hasta 3 (recta, recta por el origen o polinomio) el cálculo es exacto; si no, se hace en decimales con `np.linalg.lstsq`. Prueba subir el grado hasta que el polinomio pase por todos los datos, predecir en $x=1000$, cambiar un dato por uno atípico y usar años como $x$ (2000, 2001, ...) para ver el aviso de mal condicionamiento.

# %%
def calculadora_ajuste(x_txt, y_txt, modelo, grado, x0_txt, n_cifras):
    with bloque():                 # el texto aparece completo, de una sola vez; las gráficas van después
        try:
            res = ajuste(x_txt, y_txt, modelo, grado if modelo == "polinomio" else None, n_cifras=n_cifras)
        except ValueError as err:
            print("Revisa la entrada:", err); return
        pred, problema = None, None
        if str(x0_txt).strip():          # un problema con x0 no borra el ajuste: se dice al final
            try:
                pred = predice(res, x0_txt)
            except _Desborde as err:
                problema = f"Aviso: {err}."
            except ValueError as err:
                problema = f"Revisa «predecir en x =»: {err}"
        _informe_ajuste(res, n_cifras, pred)
        if problema:
            print(problema)
    dibuja_ajuste(res, "x", "y", "", pred, n_cifras)


_CONTROLES_4A = dict(
    x_txt=widgets.Textarea(value="0 100 200 300 400", description="x:", continuous_update=False,
                           layout=widgets.Layout(width="420px", height="60px")),
    y_txt=widgets.Textarea(value="10.0 12.1 13.8 16.2 17.9", description="y:", continuous_update=False,
                           layout=widgets.Layout(width="420px", height="60px")),
    modelo=widgets.Dropdown(options=_MODELOS, value="recta", description="modelo:"),
    grado=widgets.IntSlider(value=2, min=1, max=6, description="grado", continuous_update=False),
    x0_txt=widgets.Text(value="", description="predecir en x =", continuous_update=False,
                        style={"description_width": "initial"}),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_ajuste, **_CONTROLES_4A);

# %% [markdown]
# **Calculadora del plano.** Escribe un dato por renglón como `x, y, z`. La calculadora ajusta $z=\beta_0+\beta_1x+\beta_2y$ con las columnas $1$, $x$ y $y$ (exacto con hasta 12 datos) y dibuja el plano con los residuos.

# %%
def _tabla_plano(tabla_txt):
    """Las tres listas x, y, z (como se escribieron) de una tabla con un dato por renglón."""
    if not str(tabla_txt).strip():
        raise ValueError("escribe un dato por renglón como x, y, z, por ejemplo: 0, 0, 1.0")
    try:
        T = matriz(tabla_txt, "la tabla de datos")
    except ValueError as err:
        if str(err).startswith("escribe "):
            raise ValueError("escribe un dato por renglón como x, y, z, por ejemplo: 0, 0, 1.0") from None
        raise
    if T.cols != 3:
        raise ValueError(f"cada renglón es un dato con tres números, x, y, z (escribiste "
                         f"{_cuenta(T.cols, 'número', 'números')} por renglón), por ejemplo: 0, 0, 1.0")
    if T.rows > _MAX_DATOS:
        raise ValueError(f"la tabla tiene {T.rows} datos; esta calculadora acepta hasta {_MAX_DATOS}")
    filas, _ = _trozos(tabla_txt)
    return [f[0] for f in filas], [f[1] for f in filas], [f[2] for f in filas]


def calculadora_plano(tabla_txt, n_cifras):
    with bloque():
        try:
            res = ajuste_plano(*_tabla_plano(tabla_txt), n_cifras=n_cifras)
        except ValueError as err:
            print("Revisa la entrada:", err); return
        _aviso_coma(tabla_txt)
        _informe_ajuste(res, n_cifras)
    dibuja_ajuste(res, n_cifras=n_cifras)


_CONTROLES_4B = dict(
    tabla_txt=widgets.Textarea(value="0, 0, 1.0\n2, 0, 2.1\n0, 2, 4.9\n2, 2, 6.0\n1, 1, 3.6", description="x, y, z:",
                               continuous_update=False, layout=widgets.Layout(width="300px", height="150px")),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_plano, **_CONTROLES_4B);

# %%
# Casos de prueba de la sección 4 (resultado conocido)
_semilla_t = ajuste("0 100 200 300 400", "10.0 12.1 13.8 16.2 17.9", "recta")
_cerca = lambda a, b, tol=1e-9: abs(float(a) - float(b)) <= tol
_x14, _y14 = list(range(14)), [1, 3, 2, 4, 5, 5, 6, 8, 8, 9, 11, 12, 12, 14]       # rama decimal: más de 12 datos
_x13 = list(range(13))
_y13 = [(3 * i) % 5 for i in _x13]                                                  # puntos (x, y) no colineales
_z13 = [1 + 2 * a - b for a, b in zip(_x13, _y13)]                                  # z = 1 + 2x - y exacto
_anios, _crece = "2000 2001 2002 2003 2004 2005", "3 4 6 9 13 18"
PRUEBAS_4 = [
    ("recta por (0, 1), (1, 3), (2, 4), (3, 4) → (β0, β1) = (3/2, 1) exactos; residuos (-1/2, 1/2, 1/2, -1/2)",
     lambda: (lambda r: r["coeficientes"] == [sp.Rational(3, 2), 1] and r["exacto"]
                        and r["residuos"] == [sp.Rational(-1, 2), sp.Rational(1, 2), sp.Rational(1, 2), sp.Rational(-1, 2)])(
         ajuste("0 1 2 3", "1 3 4 4", "recta"))),
    ("mismos datos: Aᵀr = 0, AᵀA = [[4, 6], [6, 14]] y Aᵀy = (12, 23)",
     lambda: (lambda r: _es_nulo(r["A"].T * sp.Matrix(r["residuos"])) and r["ATA"] == sp.Matrix([[4, 6], [6, 14]])
                        and r["ATy"] == vector("12, 23"))(ajuste("0 1 2 3", "1 3 4 4", "recta"))),
    ("semilla del resorte (masa en gramos) → (β0, β1) = (10.02, 0.0199), SSE = 0.099 y «y = 10.02 + 0.0199 x»",
     lambda: _cerca(_semilla_t["coeficientes"][0], 10.02) and _cerca(_semilla_t["coeficientes"][1], 0.0199)
             and _semilla_t["SSE"] == sp.Rational(99, 1000) and _semilla_t["ecuacion"] == "y = 10.02 + 0.0199 x"),
    ("la recta de la semilla pasa por (x̄, ȳ) = (200, 14)",
     lambda: _semilla_t["coeficientes"][0] + 200 * _semilla_t["coeficientes"][1] == 14
             and sum(_semilla_t["x"]) / 5 == 200 and sum(_semilla_t["y"]) / 5 == 14),
    ("recta por el origen con (1, 2), (2, 4), (3, 7) → β1 = 31/14 y R2 = None",
     lambda: (lambda r: r["coeficientes"] == [sp.Rational(31, 14)] and r["R2"] is None)(
         ajuste("1 2 3", "2 4 7", "recta por el origen"))),
    ("polinomio de grado 2 con (0, 1), (1, 2), (2, 5), (3, 10) → (1, 0, 1) y SSE = 0",
     lambda: (lambda r: r["coeficientes"] == [1, 0, 1] and r["SSE"] == 0)(ajuste("0 1 2 3", "1 2 5 10", "polinomio", 2))),
    ("polinomio de grado 3 con 4 datos → SSE = 0 y el aviso de interpolación",
     lambda: (lambda r: r["SSE"] == 0 and "el polinomio pasa por todos los datos; residuo cero no significa mejor modelo"
                        in r["avisos"])(ajuste("0 1 2 3", "1 3 2 5", "polinomio", 3))),
    ("exponencial con y = 2·e^(0.5 x) exacto, x = 0, ..., 4 → a ≈ 2 y k ≈ 0.5 (tolerancia 1e-9)",
     lambda: (lambda r: _cerca(r["coeficientes"][0], 2) and _cerca(r["coeficientes"][1], 0.5))(
         ajuste([0, 1, 2, 3, 4], [2 * math.exp(0.5 * t) for t in range(5)], "exponencial"))),
    ("exponencial con un y = 0 → ValueError (se necesita y > 0)",
     lambda: "positivos" in _mensaje(lambda: ajuste("0 1 2", "1 0 2", "exponencial"))),
    ("recta con todos los x iguales → ValueError «las columnas de A son dependientes»",
     lambda: _mensaje(lambda: ajuste("2 2 2", "1 2 3", "recta")).startswith("las columnas de A son dependientes")),
    ("compara_bibliotecas con la semilla, grado 1: np.polyfit al revés = coeficientes propios (tolerancia 1e-9)",
     lambda: (lambda c: np.allclose(c["polyfit"][::-1], c["propios"], rtol=0, atol=1e-9) and c["coinciden"]
                        and not np.allclose(c["polyfit"], c["propios"], rtol=0, atol=1e-3))(
         compara_bibliotecas("0 100 200 300 400", "10.0 12.1 13.8 16.2 17.9", 1))),
    ("predice con la semilla: en x = 1000 avisa extrapolación; en x = 250 no avisa y da 14.995",
     lambda: "extrapolación" in (predice(_semilla_t, 1000)["aviso"] or "")
             and predice(_semilla_t, 250)["aviso"] is None and predice(_semilla_t, 250)["y"] == sp.Rational(14995, 1000)),
    ("datos exactamente en y = 3 - 2x → R2 = 1",
     lambda: ajuste("0 1 2 5", "3 1 -1 -7", "recta")["R2"] == 1),
    ("ajuste_plano con z = 1 + 2x - y en 4 puntos no colineales → (1, 2, -1)",
     lambda: ajuste_plano("0 1 0 1", "0 0 1 1", "1 3 0 2")["coeficientes"] == [1, 2, -1]),
    ("potencia con y = 3·x² → a ≈ 3 y b ≈ 2",
     lambda: (lambda r: _cerca(r["coeficientes"][0], 3) and _cerca(r["coeficientes"][1], 2))(
         ajuste("1 2 3 4 5", "3 12 27 48 75", "potencia"))),
    ("prueba de humo: las dos calculadoras corren con sus valores iniciales; la de ajuste, con cada modelo, con una "
     "predicción y con el deslizador en grado 4 (rama decimal, que sí usa el grado)",
     lambda: _humo(calculadora_ajuste, _CONTROLES_4A) and _humo(calculadora_plano, _CONTROLES_4B)
             and all(_humo(calculadora_ajuste, _CONTROLES_4A, modelo=m, x0_txt="5") for m in _MODELOS if m != "potencia")
             and _humo(calculadora_ajuste, _CONTROLES_4A, modelo="potencia", x_txt="1 2 3 4 5", x0_txt="5")
             and (lambda t, e: e is None and "polinomio de grado 4" in t and "x⁴" in t and "decimales" in t)(
                 *_captura(calculadora_ajuste, "0 100 200 300 400", "10.0 12.1 13.8 16.2 17.9", "polinomio", 4, "", 4))),
    ("rama decimal: polinomio de grado 4 con 6 datos y de grado 3 con 14 datos → los coeficientes de np.polyfit, al revés",
     lambda: all((lambda r, g, x, y: not r["exacto"] and np.allclose(
         [float(c) for c in r["coeficientes"]], np.polyfit(x, y, g)[::-1], rtol=1e-8, atol=1e-10))(
             ajuste(x, y, "polinomio", g), g, x, y)
         for g, x, y in ((4, [0, 1, 2, 3, 4, 5], [1, 3, 2, 5, 4, 6]), (3, _x14, _y14)))),
    ("plano en decimales (13 datos sobre z = 1 + 2x - y) → (β0, β1, β2) = (1, 2, -1) en ese orden",
     lambda: (lambda r: not r["exacto"] and np.allclose([float(c) for c in r["coeficientes"]], [1, 2, -1],
                                                       rtol=0, atol=1e-9))(ajuste_plano(_x13, _y13, _z13))),
    ("AᵀA y Aᵀy de la rama decimal (recta con 13 datos) → [[13, 78], [78, 650]] y (Σy, Σxy)",
     lambda: (lambda r: not r["exacto"] and np.allclose(r["ATA"], [[13, 78], [78, 650]])
                        and np.allclose(r["ATy"], [sum(_y13), sum(a * b for a, b in zip(_x13, _y13))]))(
         ajuste(_x13, _y13, "recta"))),
    ("RMSE y R²: en la semilla, RMSE = √(0.099/5) y R² = 1 - 0.099/39.7 exactos; en la rama decimal (14 datos), "
     "R² = 1 - SSE/SST calculado aparte",
     lambda: _semilla_t["RMSE"] == sp.sqrt(sp.Rational(99, 1000) / 5)
             and _semilla_t["R2"] == 1 - sp.Rational(99, 1000) / sp.Rational(397, 10)
             and (lambda r: _cerca(r["R2"], 1 - r["SSE"] / float(((np.array(_y14) - np.mean(_y14)) ** 2).sum()), 1e-12)
                            and _cerca(r["RMSE"], math.sqrt(r["SSE"] / 14), 1e-12))(ajuste(_x14, _y14, "polinomio", 3))),
    ("extrapolación hacia abajo: predice con la semilla en x = -100 avisa",
     lambda: "extrapolación" in (predice(_semilla_t, -100)["aviso"] or "")),
    ("años 2000 a 2005 con grado 2 (rama exacta): aviso de A mal condicionada, forma centrada en t = x - 2003, "
     "ecuación en x que reproduce el modelo (con más cifras) y AᵀA con 12015 exacto (no 12020)",
     lambda: (lambda r, t: any("mal condicionada" in a for a in r["avisos"]) and r["exacto"]
                           and "t = x - 2003" in (r["ecuacion_centrada"] or "") and r["ecuacion_reproduce"]
                           and r["cifras_ecuacion"] > 4 and "12015" in t and "12020" not in t)(
         ajuste(_anios, _crece, "polinomio", 2),
         _captura(calculadora_ajuste, _anios, _crece, "polinomio", 2, "", 4)[0])),
    ("con x en años: la exponencial avisa que ln a no cabe (sin error) y da la forma con t; la potencia aconseja cambiar "
     "la unidad de x, no restar el promedio",
     lambda: (lambda e, p: any("ln a" in a and "otro origen" in a for a in e["avisos"]) and e["ecuacion_centrada"]
                           and any("unidad" in a for a in p["avisos"]) and not any("promedio" in a for a in p["avisos"]))(
         ajuste(_anios, _crece, "exponencial"), ajuste(_anios, _crece, "potencia"))),
    ("un problema con «predecir en x =» no borra el ajuste: con «abc» y con un desborde (exponencial, x0 = 1e15)",
     lambda: (lambda a, b: "Ecuación:" in a and "Revisa «predecir en x =»" in a and "Revisa la entrada" not in a
                           and "Ecuación:" in b and "se desborda" in b)(
         _captura(calculadora_ajuste, "0 100 200 300 400", "10.0 12.1 13.8 16.2 17.9", "recta", 2, "abc", 4)[0],
         _captura(calculadora_ajuste, "0 100 200 300 400", "10.0 12.1 13.8 16.2 17.9", "exponencial", 2, "1e15", 4)[0])),
    ("compara_bibliotecas da cond(A) y cond(AᵀA), y con la semilla cond(AᵀA) ≈ cond(A)²",
     lambda: (lambda c: abs(c["cond_ATA"] / c["cond_A"] ** 2 - 1) < 1e-3 and "cond(AᵀA)" in c["texto"])(
         compara_bibliotecas("0 100 200 300 400", "10.0 12.1 13.8 16.2 17.9", 1))),
    ("plano con los puntos (x, y) casi sobre una recta: aviso, y la ecuación con las cifras que la hacen reproducir "
     "el modelo",
     lambda: (lambda r: any("casi sobre una recta" in a for a in r["avisos"]) and r["ecuacion_reproduce"]
                        and r["cifras_ecuacion"] > 4)(ajuste_plano("0 1 2 3", "0 1 2.000001 3", "1 2 3 4.1"))),
]
corre_pruebas(PRUEBAS_4, 4)

# %% [markdown]
# **Contrasta (sección 4).** Calcula a mano la recta de mínimos cuadrados $y=\beta_0+\beta_1x$ para los puntos $(0,0)$, $(1,1)$ y $(2,1)$: arma $A^TA$ y $A^T\mathbf y$ y resuelve, o usa $\beta_1=S_{xy}/S_{xx}$ y $\beta_0=\bar y-\beta_1\bar x$. Escribe tus dos coeficientes en el orden $(\beta_0,\beta_1)$, por ejemplo `mi_recta = (2, -0.5)`; puedes usar fracciones como texto (`"1/3"`) o decimales.

# %%
mi_recta = None     # dos números (β0, β1), por ejemplo (2, -0.5) o ("1/3", 2)

if mi_recta is None:     # la referencia se muestra después de tu cálculo
    print("Falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda.")
else:
    try:
        mia = vector(mi_recta, "tu recta")
    except ValueError as err:
        print(f"Revisa tu recta: {str(err).rstrip('.')}.")
    else:
        if mia.rows != 2:          # respuesta imposible: la referencia todavía no se muestra
            print(f"La recta tiene dos coeficientes, β0 y β1 (escribiste {mia.rows}): revisa tu respuesta y vuelve a "
                  "ejecutar la celda.")
        else:
            ref = ajuste("0 1 2", "0 1 1", "recta")
            b = ref["coeficientes"]
            print(f"La calculadora da: {_eq('β0', _coef_txt(b[0], 4))}, {_eq('β1', _coef_txt(b[1], 4))}  "
                  f"({ref['ecuacion']}).")
            if all(abs(float(a - c)) <= 0.01 for a, c in zip(mia, b)):
                print("Coinciden (con tolerancia de 0.01).")
            elif all(abs(float(a - c)) <= 0.01 for a, c in zip(mia, reversed(b))):
                print("NO coinciden: probablemente invertiste el orden. β0 es la ordenada al origen (el valor de y en "
                      "x = 0) y β1 es la pendiente; escribe (β0, β1).")
            else:
                print("NO coinciden: revisa AᵀA (n, Σx; Σx, Σx²) y Aᵀy (Σy, Σxy), o calcula x̄, ȳ, Sxy y Sxx.")
            print("\nEl procedimiento completo:")
            _informe_ajuste(ref, 4)
