# ID: ALG-U6-NB02
# Notebook: alg/u6_aplicaciones.ipynb · sección 6.2 análisis de mallas con leyes de Kirchhoff
# Repositorio: alg/u6_aplicaciones/02_mallas.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Análisis de mallas con leyes de Kirchhoff
#
# **Leyes de Kirchhoff.** *Ley de voltajes:* al recorrer un lazo cerrado, la suma de las subidas y las caídas de voltaje es cero. *Ley de corrientes:* la corriente que entra a un nodo es igual a la que sale.
#
# **Corrientes de malla.** A cada malla (un lazo que no encierra a otro) se le asigna una corriente $i_k$ en **sentido horario**. Con $b$ ramas y $n$ nodos hay $b-n+1$ mallas independientes. La ley de voltajes en cada malla da un sistema $Z\,\mathbf i=\mathbf V$:
#
# - $Z_{ii}$ es la suma de las resistencias de la malla $i$, y $Z_{ij}=-$(la resistencia que comparten las mallas $i$ y $j$). Por eso $Z$ es simétrica.
# - $V_i$ es la suma de las fuentes de la malla $i$: positiva si al recorrer la malla en el sentido de su corriente se pasa de $-$ a $+$ (el voltaje sube, la fuente empuja en ese sentido), negativa si se pasa de $+$ a $-$.
# - Por una rama compartida pasa la **diferencia** de las dos corrientes de malla: $i_1-i_2$ en el sentido de la malla 1.
# - Una corriente negativa no es un error: circula en sentido contrario al supuesto.
#
# **Balance de potencia.** La potencia que disipan las resistencias, $\sum R\,i^2$, es igual a la que entregan las fuentes, $\sum V\,i$. Ojo: con las corrientes que salen del mismo sistema, el balance **siempre** cierra, porque $\mathbf i^TZ\,\mathbf i=\mathbf i^T\mathbf V$ es la misma ecuación multiplicada por $\mathbf i^T$. Por eso **no detecta un signo mal puesto en una fuente**: con el signo equivocado salen otras corrientes, y el balance también cierra.
#
# **Análisis de nodos.** La alternativa es tomar como incógnitas los voltajes de los nodos respecto a tierra: con las conductancias $G=1/R$, la ley de corrientes en cada nodo da $G\,\mathbf v=\mathbf i_s$, con $G_{ii}$ = suma de las conductancias conectadas al nodo $i$, $G_{ij}=-$(la conductancia entre $i$ y $j$) e $\mathbf i_s$ las corrientes de las fuentes que entran a cada nodo.
#
# En esta sección solo hay resistencias y fuentes independientes: de voltaje en el análisis de mallas y de corriente en el de nodos. **Cómo se escriben:** en mallas, cada rama es `R, malla, otra malla` (0 si la resistencia es solo de una malla) y cada fuente es `V, malla`, con $V$ positivo si empuja en el sentido horario de esa malla; una fuente en una rama compartida se escribe dos veces, `V, 1` y `-V, 2`. En nodos, cada rama es `R, nodo, otro nodo` (0 es tierra) y cada fuente es `I, nodo`, con $I$ positiva si entra al nodo.

# %%
_MAX_MALLAS = 12                 # sistemas de hasta 12×12, exactos
_MAX_RAMAS = 30


def _ramas(ramas, que):
    """Las ramas (R, a, b) con R > 0, a de 1 a 12 y b de 0 a 12, b distinto de a. que: 'malla' o 'nodo'."""
    filas = tabla(ramas, "las ramas", "la rama", 3, _MAX_RAMAS,
                  "2, 1, 0 (R en ohms, " + ("malla y otra malla, o 0 si es solo de esa malla)" if que == "malla" else
                                            "nodo y otro nodo, o 0 si va a tierra)"))
    salida = []
    for j, (R, a, b) in enumerate(filas, 1):
        if R <= 0:
            raise ValueError(f"la rama {j} tiene R = {_ntxt(R)} Ω: una resistencia debe ser positiva (R > 0)")
        if que == "malla":
            a = _entero(a, f"en la rama {j}, el número de la malla", 1, _MAX_MALLAS)
            b = _entero(b, f"en la rama {j}, el número de la segunda malla (0 si no hay otra)", 0, _MAX_MALLAS)
            if a == b:
                raise ValueError(f"la rama {j} va de la malla {a} a la misma malla {a}: escribe dos mallas distintas (o 0 "
                                 "en la segunda si la resistencia es solo de esa malla)")
        else:
            a = _entero(a, f"en la rama {j}, el número del primer nodo (0 es tierra)", 0, _MAX_MALLAS)
            b = _entero(b, f"en la rama {j}, el número del segundo nodo (0 es tierra)", 0, _MAX_MALLAS)
            if a == b:
                raise ValueError(f"la rama {j} va del nodo {a} al mismo nodo {a}: escribe dos nodos distintos")
        salida.append((R, a, b))
    return salida


def _fuentes(fuentes, que, letra):
    """Las fuentes (valor, malla o nodo); una lista vacía o un texto vacío es «sin fuentes»."""
    if fuentes is None or (isinstance(fuentes, str) and not fuentes.strip()) or (
            isinstance(fuentes, (list, tuple)) and len(fuentes) == 0):
        return []
    filas = tabla(fuentes, "las fuentes", "la fuente", 2, _MAX_RAMAS,
                  f"10, 1 ({letra} y el número de {que})")
    return [(v, _entero(k, f"en la fuente {j}, el número de {que}", 1, _MAX_MALLAS)) for j, (v, k) in enumerate(filas, 1)]


def _textos_rama(rs, corrientes, n):
    """El texto de la corriente de cada rama, con su sentido; los decimales, con n cifras."""
    textos = []
    for j, ((R, a, b), c) in enumerate(zip(rs, corrientes), 1):
        cuenta = f"i{a}" if not b else f"i{a} - i{b}"
        lugar = f"solo en la malla {a}" if not b else f"compartida por las mallas {a} y {b}"
        valor = lambda x: _cantidad(x, "A", n)
        if c == 0:
            sentido = f"{cuenta} = 0 A: no circula corriente"
        elif c > 0:
            sentido = f"{cuenta} = {valor(c)}, en el sentido horario de la malla {a}"
        else:
            sentido = (f"{cuenta} = {valor(c)}: circulan {valor(-c)} " +
                       (f"en sentido antihorario de la malla {a}" if not b else f"en el sentido horario de la malla {b}") +
                       " (al revés de lo supuesto)")
        textos.append(f"rama {j} (R = {_ntxt(R, cif=n)} Ω, {lugar}): {sentido}")
    return textos


def mallas(ramas, fuentes):
    """Análisis de mallas, exacto. ramas: lista de (R, malla_a, malla_b), con malla_b = 0 si la resistencia es solo
    de malla_a (las mallas se numeran desde 1); fuentes: lista de (V, malla), con V positivo si empuja en el sentido
    horario de esa malla. Diccionario con Z, V, i (exactas), corrientes (la de cada rama, en el sentido horario de
    su primera malla), corrientes_rama (texto por resistencia, con su sentido), potencia_R (Σ R i²), potencia_fuentes
    (Σ V i), avisos y ramas, fuentes. Si R ≤ 0, ValueError; si Z es singular, Diagnostico (un ValueError) con la
    explicación."""
    rs, fs = _ramas(ramas, "malla"), _fuentes(fuentes, "malla", "V en volts")
    m = max([max(a, b) for _, a, b in rs] + [k for _, k in fs])
    Z, V = sp.zeros(m, m), sp.zeros(m, 1)
    for R, a, b in rs:
        Z[a - 1, a - 1] += R
        if b:
            Z[b - 1, b - 1] += R
            Z[a - 1, b - 1] -= R
            Z[b - 1, a - 1] -= R
    for v, k in fs:
        V[k - 1] += v
    vacias = [k + 1 for k in range(m) if Z[k, k] == 0]
    if vacias:
        raise ValueError(f"la malla {vacias[0]} no tiene ninguna resistencia: su renglón de Z es cero y Z es singular "
                         "(revisa la numeración de las mallas: deben ir de 1 en 1, sin saltarse ninguna)")
    if Z.det() == 0:
        raise Diagnostico("Z es singular: las ecuaciones de malla no determinan las corrientes (revisa que cada grupo de "
                         "mallas tenga al menos una resistencia que no sea compartida)")
    i = Z.LUsolve(V)
    corriente = lambda a, b: i[a - 1] - (i[b - 1] if b else 0)
    corrientes = [corriente(a, b) for _, a, b in rs]
    avisos = []
    if not fs or all(v == 0 for v, _ in fs):
        avisos.append("no hay fuentes (o todas valen 0): todas las corrientes son cero")
    return {"Z": Z, "V": V, "i": i, "corrientes": corrientes, "corrientes_rama": _textos_rama(rs, corrientes, 4),
            "potencia_R": sum(R * c ** 2 for (R, _, _), c in zip(rs, corrientes)),
            "potencia_fuentes": sum(v * i[k - 1] for v, k in fs), "avisos": avisos, "ramas": rs, "fuentes": fs}


def nodos(ramas, fuentes):
    """Análisis de nodos, exacto. ramas: lista de (R, nodo_a, nodo_b), con 0 = tierra; fuentes: lista de (I, nodo),
    con I positiva si entra al nodo. Diccionario con G (conductancias), I, v (voltajes de los nodos respecto a
    tierra, exactos), corrientes (la de cada rama, de nodo_a a nodo_b), ramas y fuentes. Si G es singular,
    ValueError."""
    rs, fs = _ramas(ramas, "nodo"), _fuentes(fuentes, "nodo", "I en amperes")
    if not any(a == 0 or b == 0 for _, a, b in rs):
        raise Diagnostico("ninguna rama llega a tierra (nodo 0): los voltajes no quedan determinados")
    m = max([max(a, b) for _, a, b in rs] + [k for _, k in fs])
    G, I = sp.zeros(m, m), sp.zeros(m, 1)
    for R, a, b in rs:
        g = 1 / R
        for k in (a, b):
            if k:
                G[k - 1, k - 1] += g
        if a and b:
            G[a - 1, b - 1] -= g
            G[b - 1, a - 1] -= g
    for c, k in fs:
        I[k - 1] += c
    sueltos = [k + 1 for k in range(m) if G[k, k] == 0]
    if sueltos:
        raise ValueError(f"el nodo {sueltos[0]} no está conectado a ninguna resistencia: su renglón de G es cero y G es "
                         "singular (revisa la numeración de los nodos)")
    if G.det() == 0:
        raise Diagnostico("G es singular: algún grupo de nodos no tiene camino a tierra por resistencias, así que sus "
                         "voltajes no quedan determinados")
    v = G.LUsolve(I)
    volt = lambda k: v[k - 1] if k else 0
    return {"G": G, "I": I, "v": v, "corrientes": [(volt(a) - volt(b)) / R for R, a, b in rs], "ramas": rs,
            "fuentes": fs}


def balance_potencia(Z, V, i):
    """Balance de potencia del análisis de mallas: potencia_R = iᵀZ i (la que disipan las resistencias) y
    potencia_fuentes = iᵀV (la que entregan las fuentes). cierra es True si son iguales."""
    Z, V, i = como_matriz(Z, "Z"), vector(V, "V"), vector(i, "i")
    if Z.rows != Z.cols:
        raise ValueError(f"Z debe ser cuadrada (una ecuación y una corriente por malla); es de {dims(Z)}")
    if V.rows != Z.rows or i.rows != Z.rows:
        raise ValueError(f"Z es de {dims(Z)}, V tiene {_cuenta(V.rows, 'número', 'números')} e i tiene {i.rows}: V e i "
                         "deben tener un número por malla")
    pR, pF = (i.T * Z * i)[0], (i.T * V)[0]
    return {"potencia_R": pR, "potencia_fuentes": pF, "cierra": pR == pF}


def _informe_mallas(res, n_cifras):
    Z, V, i = res["Z"], res["V"], res["i"]
    m = Z.rows
    print(f"{_cuenta(m, 'malla', 'mallas')}, {_cuenta(len(res['ramas']), 'rama', 'ramas')} con resistencia y "
          f"{_cuenta(len(res['fuentes']), 'fuente', 'fuentes')}. Corrientes de malla en sentido horario.")
    muestra("Z (Ω)", Z, n_cifras, tex=r"Z\ (\Omega)")
    muestra("V (V)", V, n_cifras, tex=r"\mathbf V\ (\mathrm V)")
    print("Ecuaciones de malla (ley de voltajes), Z i = V:")
    for k in range(m):
        print(f"  malla {k + 1}: {_combinacion(Z.row(k), 'i')} = {_ntxt(V[k], cif=n_cifras)}")
    print("Corrientes de malla:")
    for k in range(m):
        print("  " + _magnitud(f"i{k + 1}", i[k], "A", n_cifras) +
              ("  (negativa: circula en sentido antihorario)" if i[k] < 0 else ""))
    print("Corriente en cada rama:")
    for texto in _textos_rama(res["ramas"], res["corrientes"], n_cifras):
        print("  " + texto)
    pR, pF = res["potencia_R"], res["potencia_fuentes"]
    print(f"Balance de potencia: {_magnitud('Σ R·i²', pR, 'W', n_cifras)} y {_magnitud('Σ V·i', pF, 'W', n_cifras)}: "
          + ("cierra." if pR == pF else "¡NO cierra!"))
    for aviso in res["avisos"]:
        print("Aviso:", aviso + ".")


def _informe_nodos(res, n_cifras):
    G, I, v = res["G"], res["I"], res["v"]
    print(f"{_cuenta(G.rows, 'nodo', 'nodos')} además de tierra (0). Incógnitas: los voltajes de los nodos respecto a tierra.")
    muestra("G (S)", G, n_cifras, tex=r"G\ (\mathrm S)")
    muestra("I (A)", I, n_cifras, tex=r"\mathbf i_s\ (\mathrm A)")
    print("Voltajes de los nodos (G v = i_s):")
    for k in range(G.rows):
        print("  " + _magnitud(f"v{k + 1}", v[k], "V", n_cifras))
    print("Corriente en cada rama (de su primer nodo al segundo):")
    for j, ((R, a, b), c) in enumerate(zip(res["ramas"], res["corrientes"]), 1):
        tramo = (f"del nodo {a} " if a else "de tierra ") + (f"al nodo {b}" if b else "a tierra")
        print(f"  rama {j} (R = {_ntxt(R, cif=n_cifras)} Ω, {tramo}): "
              + _magnitud("i", c, "A", n_cifras) + ("  (circula al revés)" if c < 0 else ""))

# %% [markdown]
# **Recorrido paso a paso.** El circuito de la figura 6.3 del libro tiene dos mallas: $R_1=2\ \Omega$ solo en la malla 1, $R_3=4\ \Omega$ en la rama que comparten y $R_2=4\ \Omega$ solo en la malla 2, con la fuente $V_1=10$ V en la malla 1 (empuja en sentido horario). Con $V_2=0$:
#
# $$Z=\begin{bmatrix}6&-4\\-4&8\end{bmatrix},\qquad \mathbf i=\left(\tfrac52,\ \tfrac54\right)\ \text{A},$$
#
# lo mismo que da la reducción serie-paralelo: $2+4\parallel4=4\ \Omega$ y $10/4=2.5$ A. Con $V_2=4$ V en la malla 2 y el $+$ arriba, la fuente se opone al sentido horario de esa malla, así que entra como $-4$: $\mathbf V=(10,-4)$, $\mathbf i=(2,\ \tfrac12)$ A, y por $R_3$ pasan $i_1-i_2=1.5$ A en el sentido de la malla 1 (en la figura, hacia abajo). En la calculadora, $R_1$, $R_3$ y $R_2$ son las ramas 1, 2 y 3.

# %%
_fig63 = [(2, 1, 0), (4, 1, 2), (4, 2, 0)]          # R1, R3 (compartida) y R2 de la figura 6.3
print("(1) V1 = 10 V y V2 = 0")
_m1 = mallas(_fig63, [(10, 1)])
_informe_mallas(_m1, 4)
_paralelo = sp.Integer(4) * 4 / (4 + 4)
print(f"Serie-paralelo: R = 2 + 4∥4 = 2 + {_ntxt(_paralelo)} = {_ntxt(2 + _paralelo)} Ω  →  i1 = 10/{_ntxt(2 + _paralelo)} = "
      f"{_ntxt(10 / (2 + _paralelo))} A, y se reparte por igual entre R3 y R2: {_ntxt(5 / sp.Integer(4))} A cada una.")
print("\n" + "=" * 90 + "\n\n(2) V1 = 10 V y V2 = 4 V con el + arriba (se opone al sentido horario de la malla 2)")
_m2 = mallas(_fig63, [(10, 1), (-4, 2)])
_informe_mallas(_m2, 4)

# %% [markdown]
# **El balance de potencia no detecta un signo mal puesto.** Si en el caso (2) escribes la fuente de la malla 2 con el signo equivocado, $+4$ en lugar de $-4$, el lado derecho cambia y salen otras corrientes; aun así, el balance cierra, porque se calcula con las mismas ecuaciones. Para encontrar un signo mal puesto hay que revisar el planteamiento (el sentido en que cada fuente empuja), no el balance.

# %%
for _nombre, _fs in (("bien planteado, V = (10, -4)", [(10, 1), (-4, 2)]), ("signo equivocado, V = (10, 4)", [(10, 1), (4, 2)])):
    _r = mallas(_fig63, _fs)
    _b = balance_potencia(_r["Z"], _r["V"], _r["i"])
    print(f"{_nombre}:  i = {_tupla_u(_r['i'], 'A', 4)};  iᵀZ i = {_ntxt(_b['potencia_R'])} W,  "
          f"iᵀV = {_ntxt(_b['potencia_fuentes'])} W  →  {'cierra' if _b['cierra'] else 'no cierra'}")
print("Las dos versiones cierran el balance y dan corrientes distintas: el balance no dice cuál está bien.")

# %% [markdown]
# **Calculadora de mallas.** Escribe **una rama por renglón** como `R, malla, otra malla` (R en Ω; 0 en la otra malla si la resistencia es solo de la primera) y **una fuente por renglón** como `V, malla` (V en volts, positivo si empuja en el sentido horario de esa malla). Hasta 12 mallas. Los valores iniciales son el caso (2) del recorrido. Prueba cambiar el signo de la segunda fuente, agregar una tercera malla y escribir una resistencia de 0 Ω.

# %%
def calculadora_mallas(ramas_txt, fuentes_txt, n_cifras):
    with bloque():                 # el resultado aparece completo, de una sola vez
        try:
            res = mallas(ramas_txt, fuentes_txt)
        except ValueError as err:
            _reporta(err); return
        _aviso_coma(ramas_txt, fuentes_txt)
        _informe_mallas(res, n_cifras)


_CONTROLES_2A = dict(
    ramas_txt=widgets.Textarea(value="2, 1, 0\n4, 1, 2\n4, 2, 0", description="ramas (R en Ω, malla, malla o 0):",
                               continuous_update=False, layout=widgets.Layout(width="420px", height="110px"),
                               style={"description_width": "initial"}),
    fuentes_txt=widgets.Textarea(value="10, 1\n-4, 2", description="fuentes (V en volts, malla):",
                                 continuous_update=False, layout=widgets.Layout(width="420px", height="80px"),
                                 style={"description_width": "initial"}),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_mallas, **_CONTROLES_2A);

# %% [markdown]
# **Calculadora de nodos.** Escribe **una rama por renglón** como `R, nodo, otro nodo` (R en Ω; 0 es tierra) y **una fuente de corriente por renglón** como `I, nodo` (I en amperes, positiva si entra al nodo). Con los valores iniciales, una fuente de 3 A entra al nodo 1, que tiene 4 Ω a tierra y 2 Ω hacia el nodo 2, y el nodo 2 tiene 4 Ω a tierra.

# %%
def calculadora_nodos(ramas_txt, fuentes_txt, n_cifras):
    with bloque():
        try:
            res = nodos(ramas_txt, fuentes_txt)
        except ValueError as err:
            _reporta(err); return
        _aviso_coma(ramas_txt, fuentes_txt)
        _informe_nodos(res, n_cifras)


_CONTROLES_2B = dict(
    ramas_txt=widgets.Textarea(value="4, 1, 0\n2, 1, 2\n4, 2, 0", description="ramas (R en Ω, nodo, nodo o 0):",
                               continuous_update=False, layout=widgets.Layout(width="420px", height="110px"),
                               style={"description_width": "initial"}),
    fuentes_txt=widgets.Textarea(value="3, 1", description="fuentes (I en amperes, nodo):", continuous_update=False,
                                 layout=widgets.Layout(width="420px", height="80px"), style={"description_width": "initial"}),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_nodos, **_CONTROLES_2B);

# %%
# Casos de prueba de la sección 2 (resultado conocido)
_tres = mallas([(2, 1, 0), (3, 1, 2), (1, 2, 0), (2, 2, 3), (4, 3, 0)], [(12, 1)])
PRUEBAS_2 = [
    ("figura 6.3 con V2 = 0 → Z = [[6, -4], [-4, 8]] e i = (5/2, 5/4)",
     lambda: (lambda r: r["Z"] == sp.Matrix([[6, -4], [-4, 8]]) and r["i"] == vector("5/2, 5/4"))(
         mallas(_fig63, [(10, 1)]))),
    ("con las fuentes (10, 1) y (-4, 2) → i = (2, 1/2) y la corriente en R3 (rama 2) es 3/2",
     lambda: (lambda r: r["i"] == vector("2, 1/2") and r["corrientes"][1] == sp.Rational(3, 2)
                        and r["V"] == vector("10, -4"))(mallas("2, 1, 0\n4, 1, 2\n4, 2, 0", "10, 1\n-4, 2"))),
    ("tres mallas: ramas (2, 1, 0), (3, 1, 2), (1, 2, 0), (2, 2, 3), (4, 3, 0) y fuente (12, 1) → "
     "Z = [[5, -3, 0], [-3, 6, -2], [0, -2, 6]] e i = (192/53, 108/53, 36/53)",
     lambda: _tres["Z"] == sp.Matrix([[5, -3, 0], [-3, 6, -2], [0, -2, 6]])
             and _tres["i"] == sp.Matrix([sp.Rational(192, 53), sp.Rational(108, 53), sp.Rational(36, 53)])),
    ("Z es simétrica en las tres pruebas anteriores",
     lambda: all(r["Z"] == r["Z"].T for r in (mallas(_fig63, [(10, 1)]), mallas(_fig63, [(10, 1), (-4, 2)]), _tres))),
    ("balance de potencia: cierra con (10, 1), (-4, 2) y también con la fuente 2 de signo cambiado, con corrientes "
     "distintas",
     lambda: (lambda a, b: balance_potencia(a["Z"], a["V"], a["i"])["cierra"] and balance_potencia(b["Z"], b["V"], b["i"])["cierra"]
                           and a["potencia_R"] == a["potencia_fuentes"] and a["i"] != b["i"])(
         mallas(_fig63, [(10, 1), (-4, 2)]), mallas(_fig63, [(10, 1), (4, 2)]))),
    ("nodos([(4, 1, 0), (2, 1, 2), (4, 2, 0)], [(3, 1)]) → v = (36/5, 24/5)",
     lambda: nodos([(4, 1, 0), (2, 1, 2), (4, 2, 0)], [(3, 1)])["v"] == vector("36/5, 24/5")),
    ("una rama con R = 0 o negativa → ValueError; una malla sin resistencias → ValueError (Z singular)",
     lambda: "positiva" in _mensaje(lambda: mallas("0, 1, 0", "5, 1"))
             and "positiva" in _mensaje(lambda: mallas("2, 1, 0\n-3, 1, 2\n4, 2, 0", "5, 1"))
             and "singular" in _mensaje(lambda: mallas("2, 1, 0", "5, 1\n3, 2"))),
    ("prueba de humo: las dos calculadoras corren con sus valores iniciales, sin error y sin «Revisa la entrada»",
     lambda: _humo(calculadora_mallas, _CONTROLES_2A) and _humo(calculadora_nodos, _CONTROLES_2B)
             and _humo(calculadora_mallas, _CONTROLES_2A, fuentes_txt="10, 1\n4, 2", n_cifras=6)),
]
corre_pruebas(PRUEBAS_2, 2)

# %% [markdown]
# **Contrasta (sección 2).** Para el circuito de la figura 6.3 con $R_1=1\ \Omega$, $R_3=2\ \Omega$, $R_2=1\ \Omega$, $V_1=6$ V y $V_2=0$, escribe las dos ecuaciones de malla, resuélvelas a mano y escribe `mis_corrientes = (i1, i2)` en amperes.

# %%
mis_corrientes = None     # dos números (i1, i2) en A, por ejemplo (1.2, -0.7) o ("3/2", 1)

if mis_corrientes is None:     # la referencia se muestra después de tu cálculo
    print("Falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda.")
else:
    try:
        mias = vector(mis_corrientes, "mis_corrientes")
    except ValueError as err:
        print(f"Revisa tu respuesta: {str(err).rstrip('.')}.")
    else:
        if mias.rows != 2:          # respuesta imposible: la referencia todavía no se muestra
            print(f"El circuito tiene dos mallas, así que son dos corrientes, i1 e i2 (escribiste {mias.rows}): revisa tu "
                  "respuesta y vuelve a ejecutar la celda.")
        else:
            ref = mallas([(1, 1, 0), (2, 1, 2), (1, 2, 0)], [(6, 1)])          # los datos del enunciado
            print(f"La calculadora da: i = {_tupla_u(ref['i'], 'A', 4)}.")
            if all(abs(float(a - b)) <= 0.01 for a, b in zip(mias, ref["i"])):
                print("Coinciden (con tolerancia de 0.01 A).")
            elif all(abs(float(a - b)) <= 0.01 for a, b in zip(mias, reversed(list(ref["i"])))):
                print("NO coinciden: probablemente invertiste el orden. i1 es la corriente de la malla con la fuente.")
            else:
                print("NO coinciden: en la malla 1, Z11 = R1 + R3 y Z12 = -R3; en la malla 2, Z22 = R3 + R2. El lado "
                      "derecho es (V1, 0).")
            print("\nEl procedimiento completo:")
            _informe_mallas(ref, 4)
