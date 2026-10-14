# ID: ALG-U6-NB03
# Notebook: alg/u6_aplicaciones.ipynb · sección 6.3 mezclas y asignación de recursos
# Repositorio: alg/u6_aplicaciones/03_mezclas.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Mezclas y asignación de recursos
#
# **Mezclas.** Cada ingrediente es una **columna**: sus porcentajes de los componentes que se controlan. Las incógnitas $x_1,\dots,x_n$ son las cantidades de cada ingrediente, y hay una ecuación por el total y otra por cada componente controlado. Con dos componentes, por ejemplo,
#
# $$\begin{aligned} x_1+x_2+\dots+x_n&=\text{total}\\ p_{11}x_1+p_{12}x_2+\dots+p_{1n}x_n&=\text{total}\cdot c_1\\ p_{21}x_1+p_{22}x_2+\dots+p_{2n}x_n&=\text{total}\cdot c_2\end{aligned}\qquad\Longleftrightarrow\qquad A\mathbf x=\mathbf b ,$$
#
# con los porcentajes $p_{ij}$ y los objetivos $c_i$ escritos como fracciones (46 % es 0.46).
#
# - **La condición física** es $\mathbf x\ge\mathbf 0$: no se puede usar una cantidad negativa. Una solución con una cantidad negativa resuelve el sistema, pero la mezcla es imposible con esos ingredientes.
# - Con tantos ingredientes como ecuaciones (y composiciones independientes), la solución es **única**.
# - Con más ingredientes que ecuaciones hay **infinitas** soluciones: las cantidades dependen de variables libres, y la mezcla se escoge por otro criterio, como el costo. Escoger la más barata con $\mathbf x\ge\mathbf 0$ es un problema de *programación lineal*, que verás en un curso posterior.
# - Con tres ingredientes y dos componentes, cada ingrediente es un punto del plano de composiciones $(p_1,p_2)$, y las mezclas posibles forman el **triángulo** que tiene a los tres como vértices: una mezcla con un objetivo fuera del triángulo necesitaría una cantidad negativa.
#
# **Asignación de recursos.** Si el producto $j$ usa $a_{ij}$ horas del recurso $i$ (una máquina, una línea de ensamble) y del recurso $i$ hay $b_i$ horas, el plan que usa **todo** lo disponible cumple $A\mathbf x=\mathbf b$. Si alguna cantidad sale negativa, no se puede usar todo exactamente: hay que dejar sin usar parte de algún recurso.
#
# **Cómo se escriben.** Cada ingrediente va en un renglón con sus porcentajes, en el mismo orden de componentes que el objetivo: `46, 0` es un ingrediente con 46 % del componente 1 y 0 % del componente 2.

# %%
_MAX_INGREDIENTES, _MAX_COMPONENTES = 8, 5
_t = sp.Symbol("t")


def _intervalo(p, d):
    """Los valores de t con p + t·d ≥ 0 en todas sus componentes: (inferior, superior), con -oo u oo si no hay
    límite, o 'vacío'."""
    lo, hi = -sp.oo, sp.oo
    for pi, di in zip(p, d):
        if di > 0:
            lo = max(lo, -pi / di)
        elif di < 0:
            hi = min(hi, -pi / di)
        elif pi < 0:
            return "vacío"
    return "vacío" if lo > hi else (lo, hi)


def _analiza(A, b):
    """A x = b, exacto: estado ('única', 'infinitas' o 'sin solución'), x (si es única), solucion_general (con los
    parámetros t, o t1, t2, ...), parametros, libres (columnas libres, desde 0), intervalo (con un solo parámetro) y
    rango."""
    n = A.cols
    el = _elimina(A.row_join(b), n, con_b=True)
    base = {"x": None, "solucion_general": None, "parametros": [], "libres": [], "intervalo": None, "rango": el["r"]}
    if el["imposible"] is not None:
        return dict(base, estado="sin solución")
    R, piv = el["R"], el["pivotes"]
    xp = sp.zeros(n, 1)
    for f, c in enumerate(piv):
        xp[c] = R[f, n]
    dirs = _direcciones(R, piv, n)
    if not dirs:
        return dict(base, estado="única", x=xp)
    params = [_t] if len(dirs) == 1 else [sp.Symbol(f"t{k}") for k in range(1, len(dirs) + 1)]
    general = xp + sum((p * d for p, d in zip(params, dirs)), sp.zeros(n, 1))
    return dict(base, estado="infinitas", solucion_general=general, parametros=params,
                libres=[c for c in range(n) if c not in piv],
                intervalo=_intervalo(xp, dirs[0]) if len(dirs) == 1 else None)


def _porcentajes(valores, donde):
    for j, v in enumerate(valores, 1):
        if not 0 <= v <= 100:
            raise ValueError(f"{donde}, el porcentaje {j} vale {_ntxt(v)}: un porcentaje va de 0 a 100")
    if sum(valores) > 100:
        raise ValueError(f"{donde}, los porcentajes suman {_ntxt(sum(valores))}: no pueden pasar de 100")


def _sin_porciento(texto):
    """'20 %, 10 %' -> '20 , 10 ': el signo % se ignora (los porcentajes se escriben como números)."""
    return texto.replace("%", " ") if isinstance(texto, str) else texto


def mezcla(composiciones, objetivo, total, n_cifras=4):
    """Mezcla exacta. composiciones: un ingrediente por renglón con sus porcentajes de los componentes controlados;
    objetivo: los porcentajes deseados; total: la cantidad de mezcla. Arma A (un renglón de unos y un renglón de
    fracciones por componente) y b. Diccionario con A, b, estado, x (exacta, si es única), factible (x ≥ 0),
    avisos, solucion_general (paramétrica en t, si hay más ingredientes que ecuaciones), parametros, libres e
    intervalo (con UN parámetro libre: los t que dejan todo ≥ 0, o 'vacío'). Si dos ingredientes tienen la misma
    composición, Diagnostico (un ValueError) «... no determinan la mezcla ...». Se acepta el signo % ('20 %, 10 %').
    n_cifras: las cifras de los decimales en los avisos."""
    composiciones, objetivo = _sin_porciento(composiciones), _sin_porciento(objetivo)
    filas = tabla(composiciones, "los ingredientes", "el ingrediente", None, _MAX_INGREDIENTES,
                  "46, 0 (el porcentaje de cada componente controlado)")
    n, k = len(filas), len(filas[0])
    if k > _MAX_COMPONENTES:
        raise ValueError(f"cada ingrediente tiene {k} porcentajes; esta calculadora controla hasta {_MAX_COMPONENTES} "
                         "componentes")
    for j, f in enumerate(filas, 1):
        _porcentajes(f, f"en el ingrediente {j}")
    c = vector(objetivo, "el objetivo")
    if c.rows != k:
        raise ValueError(f"el objetivo tiene {_cuenta(c.rows, 'porcentaje', 'porcentajes')} y cada ingrediente tiene {k}: "
                         "escribe un porcentaje deseado por cada componente controlado, en el mismo orden")
    _porcentajes(list(c), "en el objetivo")
    tot = numero(total, "el total")
    if tot <= 0:
        raise ValueError(f"el total debe ser positivo (escribiste {_ntxt(tot)})")
    for i in range(n):
        for j in range(i + 1, n):
            if filas[i] == filas[j]:
                raise Diagnostico(f"los ingredientes {i + 1} y {j + 1} tienen la misma composición, así que no determinan "
                                 "la mezcla: cualquier reparto entre ellos da lo mismo (solo cuenta su suma); quita uno")
    A = sp.Matrix([[1] * n] + [[filas[j][i] / 100 for j in range(n)] for i in range(k)])
    b = sp.Matrix([tot] + [tot * c[i] / 100 for i in range(k)])
    res = _analiza(A, b)
    avisos = []
    if res["estado"] == "sin solución":
        factible = False
        avisos.append("ninguna mezcla de estos ingredientes tiene esa composición: el sistema A x = b no tiene solución")
    elif res["estado"] == "única":
        x = res["x"]
        factible = all(v >= 0 for v in x)
        for j, v in enumerate(x, 1):
            if v < 0:
                avisos.append(f"el ingrediente {j} saldría negativo: no se puede usar una cantidad negativa, así que esa "
                              "mezcla es imposible con estos ingredientes")
        cond = _condicion(A, por_renglones=True)        # renglones escalados: los porcentajes chicos no dan falsas alarmas
        if A.rows == A.cols and cond > _COND_AVISO:
            avisos.append(f"las composiciones son casi dependientes (número de condición de A, con sus renglones "
                          f"escalados, ≈ {cifras(cond, n_cifras)}): las cantidades pueden cambiar mucho con un cambio "
                          "pequeño en los porcentajes")
        if not factible and k == 1:
            lo, hi = min(f[0] for f in filas), max(f[0] for f in filas)
            avisos.append(f"el objetivo, {_ntxt(c[0])} %, está fuera del intervalo [{_ntxt(lo)}, {_ntxt(hi)}] % que "
                          "cubren los ingredientes")
    else:
        iv = res["intervalo"]
        factible = (iv != "vacío") if iv is not None else None
        if res["rango"] < min(n, k + 1):
            avisos.append(f"las composiciones son dependientes (rango de A = {res['rango']}): unos ingredientes son "
                          "mezcla de otros, así que no determinan la mezcla")
        if len(res["parametros"]) > 1:
            avisos.append("con dos o más parámetros libres, decidir si hay una mezcla con todas las cantidades no "
                          "negativas y escoger la más barata es un problema de programación lineal (curso posterior)")
        elif iv == "vacío":
            avisos.append("ningún valor del parámetro deja todas las cantidades no negativas: la mezcla es imposible")
    return dict(res, A=A, b=b, factible=factible, avisos=avisos, composiciones=filas, objetivo=c, total=tot, n=n, k=k)


def asignacion(A, b):
    """Asignación de recursos, exacta: A[i, j] = cuánto del recurso i usa una unidad del producto j; b = cuánto hay de
    cada recurso. Resuelve A x = b (usar todo lo disponible). Diccionario con A, b, estado, x, factible (x ≥ 0),
    avisos, solucion_general, parametros, libres e intervalo, como mezcla."""
    A, b = como_matriz(A, "A"), vector(b, "b")
    if A.rows != b.rows:
        raise ValueError(f"A tiene {_cuenta(A.rows, 'renglón', 'renglones')} (uno por recurso) y b tiene "
                         f"{_cuenta(b.rows, 'número', 'números')}: debe haber un dato de b por cada recurso")
    if A.rows > _MAX_EXACTO or A.cols > _MAX_EXACTO:
        raise ValueError(f"A es de {dims(A)}; esta calculadora acepta hasta {_MAX_EXACTO} recursos y {_MAX_EXACTO} productos")
    if any(e < 0 for e in A) or any(e < 0 for e in b):
        raise ValueError("los consumos de A y las disponibilidades de b no pueden ser negativos")
    res = _analiza(A, b)
    avisos, factible = [], False
    if res["estado"] == "sin solución":
        avisos.append("no hay ningún plan que use exactamente todos los recursos (el sistema no tiene solución)")
    elif res["estado"] == "única":
        factible = all(v >= 0 for v in res["x"])
        for j, v in enumerate(res["x"], 1):
            if v < 0:
                avisos.append(f"el producto {j} saldría negativo: no se puede usar todo exactamente; hay que dejar sin "
                              "usar parte de algún recurso")
    else:
        iv = res["intervalo"]
        factible = (iv != "vacío") if iv is not None else None
        if iv == "vacío":
            avisos.append("el sistema tiene infinitas soluciones, pero en todas alguna cantidad es negativa: ningún plan "
                          "con cantidades no negativas usa todo exactamente; hay que dejar sin usar parte de algún recurso")
        elif iv is None:
            avisos.append("hay infinitas soluciones; decidir cuáles tienen todas las cantidades no negativas y escoger "
                          "una por otro criterio, como la ganancia, es programación lineal (curso posterior)")
        else:
            avisos.append("hay infinitos planes que usan todo: se escoge uno por otro criterio, como la ganancia "
                          "(programación lineal, curso posterior)")
    return dict(res, A=A, b=b, factible=factible, avisos=avisos)


def _lineal_aprox(e, params, n):
    """'  ≈ 34.97 - 0.3043·t': la forma decimal (a n cifras) de una expresión lineal en los parámetros, o '' si ya es
    exacta en decimales o no tiene parte constante ni coeficientes fraccionarios."""
    const = e.subs({p: 0 for p in params})
    coefs = [e.coeff(p) for p in params]
    if all(_decimal_exacto(c, n) is not None for c in [const] + coefs):
        return ""
    partes = [] if const == 0 else [cifras(const, n)]
    for c, p in zip(coefs, params):
        if c == 0:
            continue
        signo = ("-" if c < 0 else "") if not partes else ("- " if c < 0 else "+ ")
        partes.append(f"{signo}{cifras(abs(c), n)}·{_txt(p)}")
    return "  ≈ " + " ".join(partes)


def _informe_sistema(res, letra, nombre_j, unidad, n_cifras):
    """La solución de A x = b para mezcla y asignacion: única, general o sin solución, con veredicto y comprobación."""
    A, b = res["A"], res["b"]
    if res["estado"] == "única":
        x = res["x"]
        for j, v in enumerate(x, 1):
            print("  " + _magnitud(f"{letra}{j}", v, unidad, n_cifras) + f"  ({nombre_j} {j})" +
                  ("  ← negativa" if v < 0 else ""))
        print("Comprobación: A x = b" if A * x == b else "Comprobación FALLÓ: A x ≠ b")
        print("Veredicto: es posible, todas las cantidades son no negativas (≥ 0)." if res["factible"] else
              "Veredicto: IMPOSIBLE con estos datos: el sistema tiene solución, pero alguna cantidad es negativa.")
    elif res["estado"] == "infinitas":
        params, libres = res["parametros"], res["libres"]
        print("Hay infinitas soluciones. " + ", ".join(f"{_txt(p)} = {letra}{c + 1}" for p, c in zip(params, libres)) +
              (" es libre" if len(params) == 1 else " son libres") + ":")
        for j, e in enumerate(res["solucion_general"], 1):
            print(f"  {letra}{j} = {_txt(e)}" + _lineal_aprox(e, params, n_cifras))
        iv = res["intervalo"]
        if iv == "vacío":
            print("Ningún valor del parámetro deja todas las cantidades no negativas: es imposible.")
        elif iv is not None:
            lo, hi = iv
            if lo == hi:                              # intervalo de un solo punto
                print(f"Solo un valor del parámetro deja todas las cantidades no negativas (≥ 0): t = "
                      f"{_cantidad(lo, unidad, n_cifras)}.")
            else:
                print(f"Para que todas las cantidades sean no negativas (≥ 0), t debe ir de "
                      f"{_cantidad(lo, unidad, n_cifras) if lo.is_finite else '-∞'} a "
                      f"{_cantidad(hi, unidad, n_cifras) if hi.is_finite else '∞'}.")
            for extremo in ((lo,) if lo == hi else (lo, hi)):
                if extremo.is_finite:
                    v = res["solucion_general"].subs(_t, extremo)
                    ceros = [f"{letra}{j}" for j, e in enumerate(v, 1) if e == 0]
                    print(f"  En t = {_valor(extremo, n_cifras)}: {letra} = {_tupla_u(v, unidad, n_cifras)}" +
                          (f"; {_lista(ceros)} {'vale' if len(ceros) == 1 else 'valen'} 0." if ceros else "."))
    else:
        print("El sistema A x = b no tiene solución.")
    for aviso in res["avisos"]:
        print("Aviso:", aviso + ".")


def _informe_mezcla(res, unidad, n_cifras):
    n, k, A, b = res["n"], res["k"], res["A"], res["b"]
    print(f"{_cuenta(n, 'ingrediente', 'ingredientes')}, {_cuenta(k, 'componente controlado', 'componentes controlados')}; "
          f"total {_ntxt(res['total'], cif=n_cifras)} {unidad}, con objetivo " +
          ", ".join(f"{_ntxt(v)} %" for v in res["objetivo"]) + ".")
    print("Ecuaciones (una por el total y una por cada componente):")
    print(f"  total: {_combinacion(A.row(0), 'x')} = {_ntxt(b[0])}")
    for i in range(k):
        print(f"  componente {i + 1}: {_combinacion(A.row(i + 1), 'x')} = {_ntxt(res['total'])}·{_ntxt(res['objetivo'][i] / 100)}"
              f" = {_ntxt(b[i + 1])}")
    muestra("A", A, n_cifras)
    muestra("b", b, n_cifras, tex=r"\mathbf b")
    print(f"Cantidades (en {unidad}):")
    _informe_sistema(res, "x", "ingrediente", unidad, n_cifras)


def dibuja_triangulo(composiciones, objetivos, n_cifras=4):
    """Plano de composiciones para TRES ingredientes y DOS componentes: el triángulo de las mezclas posibles (sus
    vértices son los ingredientes) y cada objetivo, dentro (círculo) o fuera (equis). Con otros tamaños no dibuja y lo
    dice."""
    filas = tabla(composiciones, "los ingredientes", "el ingrediente", None, _MAX_INGREDIENTES, "46, 0")
    if len(filas) != 3 or len(filas[0]) != 2:
        print("No se dibuja: el triángulo de mezclas es solo para tres ingredientes con dos componentes controlados.")
        return
    objs = tabla(objetivos, "los objetivos", "el objetivo", 2, 10, "20, 10")
    M = sp.Matrix([[1, 1, 1], [f[0] for f in filas], [f[1] for f in filas]])
    if M.det() == 0:
        print("No se dibuja: los tres ingredientes están alineados en el plano de composiciones (no forman un triángulo).")
        return
    V = np.array([[float(f[0]), float(f[1])] for f in filas])
    fig, ax = plt.subplots(figsize=(6.6, 7.2))
    ax.fill(V[:, 0], V[:, 1], facecolor="0.88", edgecolor="black", hatch="//", lw=1.5,
            label="mezclas posibles (todas las cantidades ≥ 0)")
    for j, (p, f) in enumerate(zip(V, filas)):
        ax.plot(*p, _MARCAS[j], color=_COLORES[j], ms=10, mec="black",
                label=f"ingrediente {j + 1}: ({_ntxt(f[0])}, {_ntxt(f[1])}) %")
    for o in objs:
        lam = M.LUsolve(sp.Matrix([1, o[0], o[1]]))
        dentro = all(v >= 0 for v in lam)
        ax.plot(float(o[0]), float(o[1]), "o" if dentro else "X", color="black" if dentro else _COLORES[1],
                ms=11 if dentro else 12, mec="black",
                label=f"objetivo ({_ntxt(o[0])}, {_ntxt(o[1])}) %: " + ("dentro" if dentro else "fuera, alguna cantidad sería negativa"))
    todos = np.vstack([V, np.array([[float(o[0]), float(o[1])] for o in objs])])
    lo, hi = todos.min(axis=0), todos.max(axis=0)
    margen = 0.12 * max(float((hi - lo).max()), 1.0)
    ax.set_xlim(lo[0] - margen, hi[0] + margen); ax.set_ylim(lo[1] - margen, hi[1] + margen)
    ax.set_aspect("equal"); ax.grid(alpha=0.4)
    ax.set_xlabel("componente 1 (%)"); ax.set_ylabel("componente 2 (%)")
    ax.set_title("Plano de composiciones: el triángulo de los tres ingredientes", fontsize=10)
    ax.legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.11), ncol=1)
    plt.tight_layout()
    plt.show()

# %% [markdown]
# **Recorrido paso a paso.** (1) Dos soluciones de un mismo ácido, al 30 % y al 70 %, para preparar 100 L al 40 %: $x_1+x_2=100$ y $0.30\,x_1+0.70\,x_2=40$ dan $(75,\ 25)$ L. Al 80 %, el sistema da $(-25,\ 125)$: imposible, porque 80 % está fuera del intervalo $[30,\,70]$ % que cubren las dos soluciones. (2) Un fertilizante de 100 kg con 20 % de nitrógeno (N) y 10 % de fósforo (P) a partir de urea (46 % N, 0 % P), fosfato diamónico (18 % N, 46 % P) y un relleno inerte (0 %, 0 %): el sistema de $3\times3$ da $\left(\tfrac{18500}{529},\ \tfrac{500}{23},\ \tfrac{22900}{529}\right)\approx(34.97,\ 21.74,\ 43.29)$ kg. En el triángulo, el objetivo $(20,10)$ queda dentro y $(30,30)$ queda fuera. (3) Un plan de dos productos: el primero usa 2 h de máquina y 1 h de ensamble, el segundo 1 h y 3 h, y hay 100 h de máquina y 150 h de ensamble: usar todo da $(30,\ 40)$ unidades.

# %%
print("(1) Dos soluciones, al 30 % y al 70 %, para 100 L al 40 %")
_informe_mezcla(mezcla("30\n70", "40", 100), "L", 4)
print("\nAl 80 %:")
_informe_mezcla(mezcla("30\n70", "80", 100), "L", 4)
print("\n" + "=" * 90 + "\n\n(2) Fertilizante: urea, fosfato diamónico y relleno")
_fert = mezcla("46, 0\n18, 46\n0, 0", "20, 10", 100)
_informe_mezcla(_fert, "kg", 4)
dibuja_triangulo("46, 0\n18, 46\n0, 0", "20, 10\n30, 30")
print("=" * 90 + "\n\n(3) Plan de producción: 2 h de máquina y 1 h de ensamble por unidad del producto 1; 1 h y 3 h del producto 2")
_plan = asignacion([[2, 1], [1, 3]], [100, 150])
muestra("A (h por unidad)", _plan["A"], tex=r"A\ (\text{h por unidad})")
muestra("b (h disponibles)", _plan["b"], tex=r"\mathbf b\ (\text{h disponibles})")
print("  máquina: " + _combinacion(_plan["A"].row(0), "x") + " = 100;   ensamble: " + _combinacion(_plan["A"].row(1), "x") + " = 150")
_informe_sistema(_plan, "x", "producto", "unidades", 4)

# %% [markdown]
# **Más ingredientes que ecuaciones.** Si al fertilizante le agregas un cuarto ingrediente con 23 % de N y 23 % de P, hay cuatro incógnitas y tres ecuaciones: la cantidad del cuarto, $t=x_4$, queda libre. Cada cantidad es una función de $t$, y pedir que las cuatro sean $\ge 0$ deja un intervalo de valores de $t$; en cada extremo, alguna cantidad vale 0. Escoger un valor de ese intervalo (por ejemplo, el más barato) ya es programación lineal.

# %%
_cuatro = mezcla("46, 0\n18, 46\n0, 0\n23, 23", "20, 10", 100)
_informe_mezcla(_cuatro, "kg", 4)

# %% [markdown]
# **Calculadora.** Escribe **un ingrediente por renglón** con sus porcentajes de los componentes controlados (en el mismo orden que el objetivo), el objetivo en porcentaje, el total y su unidad. Con tres ingredientes y dos componentes, la calculadora dibuja el triángulo de mezclas posibles. Prueba un objetivo fuera del triángulo (`30, 30`), un cuarto ingrediente (`23, 23`) y dos ingredientes iguales.

# %%
def calculadora_mezcla(ingredientes_txt, objetivo_txt, total_txt, unidad, n_cifras):
    with bloque():                 # el texto aparece completo, de una sola vez; el dibujo va después
        try:
            u = nombre_unidad(unidad, "cantidad", "kg o L")
            res = mezcla(ingredientes_txt, objetivo_txt, total_txt, n_cifras=n_cifras)
        except ValueError as err:
            _reporta(err); return
        _aviso_coma(ingredientes_txt, objetivo_txt)
        _informe_mezcla(res, u, n_cifras)
    if res["n"] == 3 and res["k"] == 2:
        dibuja_triangulo(res["composiciones"], [list(res["objetivo"])], n_cifras)


_CONTROLES_3 = dict(
    ingredientes_txt=widgets.Textarea(value="46, 0\n18, 46\n0, 0", description="ingredientes (%):",
                                      continuous_update=False, layout=widgets.Layout(width="320px", height="110px"),
                                      style={"description_width": "initial"}),
    objetivo_txt=widgets.Text(value="20, 10", description="objetivo (%):", continuous_update=False,
                              style={"description_width": "initial"}),
    total_txt=widgets.Text(value="100", description="total:", continuous_update=False),
    unidad=widgets.Text(value="kg", description="unidad:", continuous_update=False),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_mezcla, **_CONTROLES_3);

# %%
# Casos de prueba de la sección 3 (resultado conocido)
_iv = _cuatro["intervalo"]
PRUEBAS_3 = [
    ("dos soluciones, 30 % y 70 %, para 100 al 40 % → (75, 25) y factible; para 50 al 40 % → (37.5, 12.5)",
     lambda: (lambda r: r["x"] == vector("75, 25") and r["factible"] and not r["avisos"])(mezcla([(30,), (70,)], (40,), 100))
             and mezcla("30\n70", "40", 50)["x"] == vector("37.5, 12.5")),
    ("objetivo 80 % → (-25, 125), no factible y con aviso",
     lambda: (lambda r: r["x"] == vector("-25, 125") and r["factible"] is False
                        and any("ingrediente 1 saldría negativo" in a for a in r["avisos"])
                        and any("fuera del intervalo" in a for a in r["avisos"]))(mezcla("30\n70", "80", 100))),
    ("fertilizante → (18500/529, 500/23, 22900/529) exactas y factible (también si se escriben los signos %)",
     lambda: _fert["x"] == sp.Matrix([sp.Rational(18500, 529), sp.Rational(500, 23), sp.Rational(22900, 529)])
             and _fert["factible"]
             and mezcla("46 %, 0 %\n18 %, 46 %\n0, 0", "20 %, 10 %", 100)["x"] == _fert["x"]),
    ("cuatro ingredientes (con 23 %, 23 %) → solución general con un parámetro t = x4 e intervalo [0, 1000/23]; en "
     "los extremos alguna cantidad es 0 y fuera de él alguna es negativa",
     lambda: _cuatro["estado"] == "infinitas" and _cuatro["parametros"] == [_t] and _cuatro["libres"] == [3]
             and _iv == (0, sp.Rational(1000, 23)) and _cuatro["factible"]
             and all(any(e == 0 for e in _cuatro["solucion_general"].subs(_t, v)) for v in _iv)
             and all(any(e < 0 for e in _cuatro["solucion_general"].subs(_t, v)) for v in (-1, _iv[1] + 1))
             and _cuatro["A"] * _cuatro["solucion_general"].subs(_t, 7) == _cuatro["b"]
             and mezcla("0\n100\n50", "0", 10)["intervalo"] == (0, 0)
             and "Solo un valor" in _captura(calculadora_mezcla, "0\n100\n50", "0", "10", "kg", 4)[0]),
    ("asignacion([[2, 1], [1, 3]], [100, 150]) → (30, 40); con A = [[1, 2], [3, 1]] y b = (50, 60) → (14, 18) (no la "
     "de Aᵀ); [[1, 1, 1], [1, 2, 3]] con (10, 40): infinitas soluciones, ninguna no negativa",
     lambda: (lambda r: r["x"] == vector("30, 40") and r["factible"])(asignacion([[2, 1], [1, 3]], [100, 150]))
             and asignacion([[1, 2], [3, 1]], [50, 60])["x"] == vector("14, 18")
             and (lambda r: r["estado"] == "infinitas" and r["factible"] is False and r["intervalo"] == "vacío"
                            and any("ningún plan con cantidades no negativas" in a for a in r["avisos"]))(
                 asignacion([[1, 1, 1], [1, 2, 3]], [10, 40]))),
    ("dos ingredientes con la misma composición y nada más → ValueError «no determinan la mezcla» (sin traza)",
     lambda: "no determinan la mezcla" in _mensaje(lambda: mezcla("30\n30", "40", 100))
             and "no determinan la mezcla" in _captura(calculadora_mezcla, "30\n30", "30", "100", "L", 4)[0]),
    ("prueba de humo: la calculadora corre con sus valores iniciales, con dos ingredientes y con cuatro, sin error y "
     "sin «Revisa la entrada»",
     lambda: _humo(calculadora_mezcla, _CONTROLES_3)
             and _humo(calculadora_mezcla, _CONTROLES_3, ingredientes_txt="30\n70", objetivo_txt="80", unidad="L")
             and _humo(calculadora_mezcla, _CONTROLES_3, ingredientes_txt="46, 0\n18, 46\n0, 0\n23, 23", n_cifras=6)),
]
corre_pruebas(PRUEBAS_3, 3)

# %% [markdown]
# **Contrasta (sección 3).** Dos aleaciones tienen 90 % y 60 % de cobre. ¿Cuántos kg de cada una necesitas para 30 kg al 70 %? Plantea las dos ecuaciones (total y cobre), resuélvelas a mano y escribe `mis_kg = (a, b)`: primero la de 90 %.

# %%
mis_kg = None     # dos números (a, b) en kg, por ejemplo (12.5, 17.5)

if mis_kg is None:     # la referencia se muestra después de tu cálculo
    print("Falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda.")
else:
    try:
        mias = vector(mis_kg, "mis_kg")
    except ValueError as err:
        print(f"Revisa tu respuesta: {str(err).rstrip('.')}.")
    else:
        if mias.rows != 2:          # respuestas imposibles: la referencia todavía no se muestra
            print(f"Son dos aleaciones, así que son dos cantidades (escribiste {mias.rows}): revisa tu respuesta y vuelve "
                  "a ejecutar la celda.")
        elif any(v < 0 for v in mias):
            print("Una cantidad de aleación no puede ser negativa: revisa tu cálculo y vuelve a ejecutar la celda.")
        elif abs(float(sum(mias)) - 30) > 0.01:
            print(f"Tus cantidades suman {_ntxt(sum(mias))} kg y deben sumar 30 kg (es la primera ecuación): revisa y "
                  "vuelve a ejecutar la celda.")
        else:
            ref = mezcla("90\n60", "70", 30)          # los datos del enunciado
            print(f"La calculadora da: (a, b) = {_tupla_u(ref['x'], 'kg', 4)}.")
            if all(abs(float(a - b)) <= 0.01 for a, b in zip(mias, ref["x"])):
                print("Coinciden (con tolerancia de 0.01 kg).")
            elif all(abs(float(a - b)) <= 0.01 for a, b in zip(mias, reversed(list(ref["x"])))):
                print("NO coinciden: probablemente invertiste el orden; a es la aleación de 90 %.")
            else:
                print("NO coinciden: la ecuación del cobre es 0.90 a + 0.60 b = 0.70·30 = 21.")
            print("\nEl procedimiento completo:")
            _informe_mezcla(ref, "kg", 4)
