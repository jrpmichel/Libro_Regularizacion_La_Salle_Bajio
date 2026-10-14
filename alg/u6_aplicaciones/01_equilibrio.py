# ID: ALG-U6-NB01
# Notebook: alg/u6_aplicaciones.ipynb · sección 6.1 equilibrio de fuerzas en 3D, momento y trabajo
# Repositorio: alg/u6_aplicaciones/01_equilibrio.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Equilibrio de fuerzas en 3D, momento y trabajo
#
# **Fuerza a lo largo de un cable.** Un cable que va del nudo $P$ al anclaje $Q$ jala al nudo hacia $Q$. Si su tensión es $T$, la fuerza que hace sobre el nudo es
#
# $$\mathbf F=T\,\mathbf u,\qquad \mathbf u=\frac{Q-P}{\|Q-P\|}\quad(\text{el vector unitario de }P\text{ hacia }Q).$$
#
# **Equilibrio.** Una partícula (un nudo donde se juntan cables) está en equilibrio si $\sum\mathbf F=\mathbf 0$: en 3D son tres ecuaciones, una por componente. Un cuerpo rígido, que además puede girar, necesita $\sum\mathbf F=\mathbf 0$ y $\sum\mathbf M=\mathbf 0$: seis ecuaciones.
#
# **Tres cables y una carga.** Si en el nudo actúa la carga $\mathbf W$ y los vectores unitarios de los cables son las columnas de $U=[\,\mathbf u_1\ \mathbf u_2\ \mathbf u_3\,]$,
#
# $$T_1\mathbf u_1+T_2\mathbf u_2+T_3\mathbf u_3+\mathbf W=\mathbf 0\qquad\Longleftrightarrow\qquad U\mathbf T=-\mathbf W .$$
#
# - Un cable solo jala. Si alguna tensión sale negativa, ese cable tendría que empujar: **con cables, esa geometría no sostiene la carga**, aunque el sistema tenga solución. Si en lugar de cables son barras o patas (un trípode, por ejemplo), que sí pueden empujar, una tensión negativa es una **compresión** y la geometría sí sostiene la carga.
# - Si los tres cables están en un mismo plano, $U$ es singular ($\det U=0$): no hay solución única, y una carga fuera de ese plano no se puede equilibrar.
#
# **Momento.** El momento respecto a un punto $O$ de una fuerza $\mathbf F$ aplicada en $P$ es
#
# $$\mathbf M_O=\mathbf r\times\mathbf F,\qquad \mathbf r=P-O\quad(\text{del punto }O\text{ al punto de aplicación}).$$
#
# El orden importa: $\mathbf F\times\mathbf r=-\mathbf M_O$. $\mathbf M_O$ es perpendicular a $\mathbf r$ y a $\mathbf F$ (es el eje alrededor del cual la fuerza tiende a hacer girar), y su norma es la fuerza por el brazo de palanca: $\|\mathbf M_O\|=\|\mathbf F\|\,d$, con $d$ la distancia de $O$ a la línea de acción de $\mathbf F$.
#
# **Trabajo.** Una fuerza constante $\mathbf F$ a lo largo de un desplazamiento recto $\mathbf d$ hace el trabajo $W=\mathbf F\cdot\mathbf d=\|\mathbf F\|\,\|\mathbf d\|\cos\theta$ (aquí $W$ es el trabajo, no la carga): en J si $\mathbf F$ está en N y $\mathbf d$ en m. Solo cuenta la componente de la fuerza a lo largo del desplazamiento; una fuerza perpendicular al desplazamiento no hace trabajo.

# %%
def _vector3(p, nombre):
    """Un vector de R³ (tres componentes x, y, z) como columna de racionales exactos."""
    v = vector(p, nombre)
    if v.rows != 3:
        raise ValueError(f"{nombre} debe tener tres componentes, x, y, z (escribiste {_cuenta(v.rows, 'componente', 'componentes')})")
    return v


def vector_unitario(P, Q):
    """u = (Q - P)/‖Q - P‖, exacto (con raíces si hace falta): la dirección de P hacia Q, de longitud 1. Si P = Q,
    ValueError «P y Q son el mismo punto»."""
    P, Q = vector(P, "P"), vector(Q, "Q")
    if P.rows != Q.rows:
        raise ValueError(f"P tiene {_cuenta(P.rows, 'componente', 'componentes')} y Q tiene {Q.rows}: deben tener el "
                         "mismo número de componentes")
    d = Q - P
    if _es_nulo(d):
        raise ValueError("P y Q son el mismo punto: no definen ninguna dirección (‖Q - P‖ = 0)")
    return d / sp.sqrt(_punto(d, d))


def _anclajes(anclajes):
    """Los anclajes como lista de columnas exactas de R³ (uno por renglón)."""
    filas = tabla(anclajes, "los anclajes", "el anclaje", 3, 12, "2, 2, 1 (las coordenadas x, y, z del anclaje)")
    return [sp.Matrix(f) for f in filas]


def equilibrio_particula(anclajes, carga, nudo=(0, 0, 0), n_cifras=4):
    """Equilibrio de un nudo sostenido por TRES cables: resuelve U T = -W, exacto. anclajes: los tres puntos Q1, Q2,
    Q3 (uno por renglón); carga: el vector W que actúa en el nudo; nudo: el punto P donde se juntan los cables.
    Diccionario con U (vectores unitarios como columnas, exactos), tensiones (lista exacta), todas_positivas (True si
    las tres son > 0), sostiene (True si ninguna es negativa), avisos (por ejemplo, «el cable 3 tendría que empujar»),
    comprobacion (T1·u1 + T2·u2 + T3·u3 + W, que debe ser el vector cero), equilibrio (True si lo es), det_U,
    longitudes, direcciones (Qj - P), anclajes, W, nudo y cond (número de condición de U). Si los tres cables están en
    un mismo plano, Diagnostico (un ValueError) «los cables están en un mismo plano: con esta geometría no hay solución
    única». n_cifras: las cifras de los decimales en los avisos."""
    Q = _anclajes(anclajes)
    if len(Q) != 3:
        raise ValueError(f"escribiste {_cuenta(len(Q), 'anclaje', 'anclajes')}: este planteamiento usa exactamente tres "
                         "cables, porque en 3D hay tres ecuaciones (una por componente) y U debe ser cuadrada")
    W = _vector3(carga, "la carga")
    N = _vector3(nudo, "el nudo")
    D = [q - N for q in Q]
    for j, d in enumerate(D, 1):
        if _es_nulo(d):
            raise ValueError(f"el anclaje {j} coincide con el nudo: ese cable no tiene dirección")
    V = sp.Matrix.hstack(*D)                   # direcciones sin normalizar: entradas racionales
    det_V = V.det()
    if det_V == 0:
        raise Diagnostico("los cables están en un mismo plano: con esta geometría no hay solución única "
                          "(U es singular, det U = 0)")
    L = [sp.sqrt(_punto(d, d)) for d in D]
    U = sp.Matrix.hstack(*[d / l for d, l in zip(D, L)])
    q = V.LUsolve(-W)                          # V q = -W con q_j = T_j/L_j: se resuelve con fracciones, sin raíces
    T = [q[j] * L[j] for j in range(3)]
    suma = _simplifica(W + sum((T[j] * U[:, j] for j in range(3)), sp.zeros(3, 1)))
    avisos = []
    for j, t in enumerate(T, 1):
        if t < 0:
            avisos.append(f"el cable {j} tendría que empujar (T{j} < 0): un cable solo jala, así que con cables esta "
                          "geometría no sostiene la carga; si son barras o patas (pueden empujar), una tensión negativa "
                          "es compresión y la geometría sí sostiene la carga")
        elif t == 0:
            avisos.append(f"el cable {j} queda flojo (T{j} = 0): con esta carga no trabaja")
    cond = _condicion(U)
    if cond > _COND_AVISO:
        avisos.append(f"los cables están casi en un mismo plano (número de condición de U ≈ {cifras(cond, n_cifras)}): "
                      "las tensiones pueden salir muy grandes y cambian mucho si un anclaje se mueve un poco")
    return {"U": U, "tensiones": T, "todas_positivas": all(t > 0 for t in T), "sostiene": all(t >= 0 for t in T),
            "avisos": avisos, "comprobacion": suma, "equilibrio": _es_nulo(suma), "det_U": det_V / (L[0] * L[1] * L[2]),
            "longitudes": L, "direcciones": D, "anclajes": Q, "W": W, "nudo": N, "cond": cond}


def momento(F, P, O=(0, 0, 0)):
    """Momento de la fuerza F aplicada en el punto P respecto al punto O: M = r × F con r = P - O, exacto."""
    F, P, O = _vector3(F, "F"), _vector3(P, "P"), _vector3(O, "O")
    return (P - O).cross(F)


def trabajo(F, d):
    """Trabajo de una fuerza constante F en el desplazamiento recto d: W = F·d, exacto."""
    F, d = vector(F, "F"), vector(d, "d")
    if F.rows != d.rows:
        raise ValueError(f"F tiene {_cuenta(F.rows, 'componente', 'componentes')} y d tiene {d.rows}: deben tener el "
                         "mismo número de componentes")
    return _punto(F, d)


def _par_txt(e):
    """Un número dentro de una cuenta: 0.3, 50, (-200), (1/3)."""
    t = _ntxt(e)
    return f"({t})" if (e < 0 or "/" in t or "√" in t) else t


def _informe_momento(F, P, O, n_cifras, uF="N", uL="m"):
    F, P, O = _vector3(F, "F"), _vector3(P, "P"), _vector3(O, "O")
    r = P - O
    M = r.cross(F)
    (rx, ry, rz), (Fx, Fy, Fz) = list(r), list(F)
    p = _par_txt
    print(f"r = P - O = {_tupla_u(r, uL)},  F = {_tupla_u(F, uF)}")
    print("M = r × F = (ry·Fz - rz·Fy,  rz·Fx - rx·Fz,  rx·Fy - ry·Fx)")
    print(f"  = ({p(ry)}·{p(Fz)} - {p(rz)}·{p(Fy)},  {p(rz)}·{p(Fx)} - {p(rx)}·{p(Fz)},  "
          f"{p(rx)}·{p(Fy)} - {p(ry)}·{p(Fx)})")
    print(f"  = {_tupla_u(M, uF + '·' + uL, n_cifras)}")
    perpendicular = _punto(M, r) == 0 and _punto(M, F) == 0
    print(_magnitud("‖M‖", sp.sqrt(_punto(M, M)), uF + "·" + uL, n_cifras) +
          ("" if _es_nulo(M) else (";  M·r = 0 y M·F = 0: M es perpendicular a r y a F" if perpendicular else
                                   ";  ¡M no es perpendicular a r y a F!")))
    return M


def _informe_equilibrio(res, unidad, n_cifras):
    N, W, U = res["nudo"], res["W"], res["U"]
    print(f"Nudo P = {_tupla_u(N, 'm')};  carga W = {_tupla_u(W, unidad)}")
    for j in range(3):
        d, l = res["direcciones"][j], res["longitudes"][j]
        print(f"Cable {j + 1}: Q{j + 1} - P = {_tupla_u(d, 'm')},  {_magnitud(f'‖Q{j + 1} - P‖', l, 'm', n_cifras)},  "
              f"u{j + 1} = {_tupla(U[:, j], n_cifras)}")
    muestra("U", U, n_cifras, tex=r"U=[\,\mathbf u_1\ \mathbf u_2\ \mathbf u_3\,]")
    print(_magnitud("det U", res["det_U"], "", n_cifras) + ";  es distinto de 0: los cables no están en un mismo plano")
    print(f"Equilibrio del nudo: T1·u1 + T2·u2 + T3·u3 + W = 0, es decir, U T = -W = {_tupla_u(-W, unidad)}")
    if all(_es_cero(e) for e in U.T * U - sp.eye(3)):
        print("Aquí UᵀU = I (los tres cables son perpendiculares entre sí), así que T = -UᵀW sin resolver nada.")
    print("Tensiones:")
    for j, t in enumerate(res["tensiones"], 1):
        print("  " + _magnitud(f"T{j}", t, unidad, n_cifras))
    if res["todas_positivas"]:
        print("Veredicto: los tres cables jalan (las tres tensiones son positivas): esta geometría sostiene la carga.")
    elif res["sostiene"]:
        print("Veredicto: esta geometría sostiene la carga, pero no todos los cables trabajan.")
    else:
        print("Veredicto: con cables, esta geometría NO sostiene la carga (un cable no puede empujar). Si son barras o "
              "patas, que sí pueden empujar, las tensiones negativas son compresiones y la geometría sí la sostiene.")
    for aviso in res["avisos"]:
        print("Aviso:", aviso + ".")
    print("Comprobación: T1·u1 + T2·u2 + T3·u3 + W = " + _tupla_u(res["comprobacion"], unidad) +
          (" (la suma de fuerzas es cero)" if res["equilibrio"] else " (¡no es cero!)"))


def dibuja_cables(anclajes, carga, tensiones, nudo=(0, 0, 0), unidad="N", n_cifras=4):
    """Dibujo en 3D del nudo, los cables (un estilo de línea y un marcador por cable), la dirección unitaria de cada
    cable (flecha corta, no a escala), su tensión y la carga (flecha negra gruesa). tensiones puede ser None."""
    Q = _anclajes(anclajes)
    W, N = _vector3(carga, "la carga"), _vector3(nudo, "el nudo")
    Nf, Wf = _flotante(N).ravel(), _flotante(W).ravel()
    Qf = [_flotante(q).ravel() for q in Q]
    largo = max([float(np.linalg.norm(q - Nf)) for q in Qf] + [1e-12])
    fig = plt.figure(figsize=(7.4, 8.0))
    ax = fig.add_subplot(projection="3d")
    fig.subplots_adjust(bottom=0.2)
    puntos = [Nf] + Qf
    for j, q in enumerate(Qf):
        k = j % len(_COLORES)
        texto_T = ""
        if tensiones is not None:
            t = tensiones[j]
            texto_T = f"T{j + 1} = {_dec_o_cifras(t, n_cifras)} {unidad}" + ("  (empujaría)" if float(t) < 0 else "")
        ax.plot(*np.array([Nf, q]).T, linestyle=_ESTILOS[k], color=_COLORES[k], lw=2.2,
                label=f"cable {j + 1}" + (f": {texto_T}" if texto_T else ""))
        ax.scatter(*q, marker=_MARCAS[k], s=60, color=_COLORES[k], depthshade=False)
        d = q - Nf
        if np.linalg.norm(d) > 1e-12 * largo:
            u = d / np.linalg.norm(d)
            ax.quiver(*Nf, *(0.3 * largo * u), color=_COLORES[k], linewidth=2.8, arrow_length_ratio=0.3)
        ax.text(*q, f" Q{j + 1}", fontsize=9)
        if texto_T:                                   # la tensión, a media altura del cable (lejos de los ejes)
            ax.text(*(Nf + 0.6 * d), f" T{j + 1}", fontsize=9, color="black",
                    bbox=dict(boxstyle="round,pad=0.1", facecolor="white", edgecolor="none", alpha=0.75))
    nW = float(np.linalg.norm(Wf))
    if nW > 0:
        punta = Nf + 0.45 * largo * Wf / nW
        ax.quiver(*Nf, *(punta - Nf), color="black", linewidth=3.2, arrow_length_ratio=0.22)
        ax.plot([], [], "-", color="black", lw=3.2, label=f"carga W: ‖W‖ = {_dec_o_cifras(sp.sqrt(_punto(W, W)), n_cifras)} {unidad}")
        ax.text(*punta, f"  W", fontsize=10)
        puntos.append(punta)
    ax.scatter(*Nf, marker="o", s=70, color="black", depthshade=False)
    ax.text(*Nf, "  P (nudo)", fontsize=9)
    ax.plot([], [], "-", color="0.4", lw=2.8, label="flechas cortas: dirección unitaria de cada cable (no a escala)")
    Pm = np.array(puntos)
    centro = (Pm.max(axis=0) + Pm.min(axis=0)) / 2
    radio = max(float((Pm.max(axis=0) - Pm.min(axis=0)).max()) / 2, 1e-9) * 1.15
    ax.set_xlim(centro[0] - radio, centro[0] + radio)
    ax.set_ylim(centro[1] - radio, centro[1] + radio)
    ax.set_zlim(centro[2] - radio, centro[2] + radio)
    ax.set_box_aspect((1, 1, 1))                          # misma escala en los tres ejes
    ax.set_xlabel("x (m)"); ax.set_ylabel("y (m)"); ax.set_zlabel("z (m)")
    ax.view_init(elev=20, azim=-60)
    ax.set_title("Nudo en equilibrio: cables, direcciones unitarias y carga", fontsize=10)
    fig.legend(*ax.get_legend_handles_labels(), loc="lower center", ncol=2, fontsize=8)
    plt.show()

# %% [markdown]
# **Recorrido paso a paso.** La lámpara del libro: cuelga de un nudo en el origen, sostenido por tres cables con anclajes en $(2,2,1)$, $(-2,1,2)$ y $(1,-2,2)$ m, y pesa 90 N, así que la carga es $\mathbf W=(0,0,-90)$ N. Cada cable mide 3 m, los vectores unitarios son $\tfrac13(2,2,1)$, $\tfrac13(-2,1,2)$ y $\tfrac13(1,-2,2)$, y resolver $U\mathbf T=-\mathbf W$ da $\mathbf T=(30,\,60,\,60)$ N: los tres cables jalan. Después, el momento de $\mathbf F=(0,0,-200)$ N aplicada en $P=(0.3,\,0.4,\,0)$ m respecto al origen: $\mathbf M=\mathbf r\times\mathbf F=(-80,\,60,\,0)$ N·m, de norma $100=200\cdot0.5$ (la fuerza es vertical y el brazo de palanca es la distancia horizontal de $O$ a $P$, 0.5 m); con el orden invertido sale $\mathbf F\times\mathbf r=(80,\,-60,\,0)$. Al final, el trabajo de $\mathbf F=(40,30,0)$ N: en $\mathbf d=(5,0,0)$ m hace 200 J, y en $\mathbf d=(3,4,0)$ m, del mismo largo, 240 J, porque ahí la fuerza y el desplazamiento forman un ángulo menor.

# %%
_lampara = equilibrio_particula("2, 2, 1; -2, 1, 2; 1, -2, 2", "0, 0, -90")
_informe_equilibrio(_lampara, "N", 4)
dibuja_cables(_lampara["anclajes"], _lampara["W"], _lampara["tensiones"], unidad="N")
print("\n" + "=" * 90 + "\n\nMomento de F = (0, 0, -200) N aplicada en P = (0.3, 0.4, 0) m, respecto al origen:")
_M = _informe_momento("0, 0, -200", "0.3, 0.4, 0", "0, 0, 0", 4)
_brazo = sp.sqrt(sp.Rational(3, 10) ** 2 + sp.Rational(2, 5) ** 2)
print(f"Brazo de palanca: F es vertical, así que es la distancia horizontal de O a P, √(0.3² + 0.4²) = "
      f"{_ntxt(_brazo)} m;  ‖F‖·d = 200·{_ntxt(_brazo)} = {_ntxt(200 * _brazo)} N·m = ‖M‖")
print(f"Con el orden invertido: F × r = {_tupla_u(vector('0, 0, -200').cross(vector('0.3, 0.4, 0')), 'N·m')} = -M")
print("\n" + "=" * 90 + "\n\nTrabajo de F = (40, 30, 0) N:")
for _d in ("5, 0, 0", "3, 4, 0"):
    _Wt = trabajo("40, 30, 0", _d)
    print(f"  en d = {_tupla_u(vector(_d), 'm')}: W = F·d = {_txt(_Wt)} J;  cos θ = W/(‖F‖·‖d‖) = {_txt(_Wt)}/(50·5) = "
          f"{_ntxt(_Wt / 250)}")
print(f"  con d = (0, 0, 2) m, perpendicular a F: W = {_txt(trabajo('40, 30, 0', '0, 0, 2'))} J")

# %% [markdown]
# **Calculadora.** Escribe los **tres anclajes**, uno por renglón, como `x, y, z` en metros (el nudo está en el origen), la carga $\mathbf W$ que actúa en el nudo y la unidad de fuerza. La calculadora arma $U$, resuelve $U\mathbf T=-\mathbf W$, dice si la geometría sostiene la carga y dibuja los cables. Si en lugar de cables son barras o patas (pueden empujar), una tensión negativa es una compresión y la geometría sí sostiene la carga. Prueba mover un anclaje en horizontal hasta que, vistos desde arriba, dos cables queden alineados (por ejemplo, el anclaje 2 en `-2, -2, 2`: sale $\mathbf T=(90,\ 60\sqrt3,\ 0)$ N y el cable 3 queda flojo), poner los tres anclajes a la misma altura que el nudo (los cables quedan en un plano) y los anclajes `1, 0, 1`, `-1, 0, 1`, `0, 1, 1` con la carga `0, 2, -10`.

# %%
def calculadora_equilibrio(anclajes_txt, carga_txt, unidad, n_cifras):
    with bloque():                 # el texto aparece completo, de una sola vez; el dibujo va después
        try:
            u = nombre_unidad(unidad, "fuerza", "N o kN")
            res = equilibrio_particula(anclajes_txt, carga_txt, n_cifras=n_cifras)
        except Diagnostico as err:            # cables en un plano: se dice y se dibuja la geometría, sin tensiones
            _reporta(err)
            res = None
        except ValueError as err:
            _reporta(err); return
        _aviso_coma(anclajes_txt, carga_txt)
        if res is not None:
            _informe_equilibrio(res, u, n_cifras)
    if res is None:
        dibuja_cables(anclajes_txt, carga_txt, None, unidad=u, n_cifras=n_cifras)
    else:
        dibuja_cables(res["anclajes"], res["W"], res["tensiones"], unidad=u, n_cifras=n_cifras)


_CONTROLES_1 = dict(
    anclajes_txt=widgets.Textarea(value="2, 2, 1\n-2, 1, 2\n1, -2, 2", description="anclajes (m):",
                                  continuous_update=False, layout=widgets.Layout(width="320px", height="100px"),
                                  style={"description_width": "initial"}),
    carga_txt=widgets.Text(value="0, 0, -90", description="carga W =", continuous_update=False),
    unidad=widgets.Text(value="N", description="unidad de fuerza:", continuous_update=False,
                        style={"description_width": "initial"}),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False))
widgets.interact(calculadora_equilibrio, **_CONTROLES_1);

# %%
# Casos de prueba de la sección 1 (resultado conocido)
_r2 = sp.sqrt(2)
PRUEBAS_1 = [
    ("vector_unitario((0, 0, 0), (2, 2, 1)) = (2/3, 2/3, 1/3)",
     lambda: tuple(vector_unitario((0, 0, 0), (2, 2, 1))) == (sp.Rational(2, 3), sp.Rational(2, 3), sp.Rational(1, 3))),
    ("lámpara: anclajes (2, 2, 1), (-2, 1, 2), (1, -2, 2) y carga (0, 0, -90) → T = (30, 60, 60) exactas, todas "
     "positivas y suma de fuerzas cero",
     lambda: (lambda r: r["tensiones"] == [30, 60, 60] and r["todas_positivas"] and r["equilibrio"]
                        and r["comprobacion"] == sp.zeros(3, 1) and not r["avisos"])(
         equilibrio_particula("2, 2, 1; -2, 1, 2; 1, -2, 2", "0, 0, -90"))),
    ("anclajes (1, 0, 1), (-1, 0, 1), (0, 1, 1) con carga (0, 2, -10) → T = (6√2, 6√2, -2√2) y el aviso del cable 3",
     lambda: (lambda r: [sp.simplify(t - e) for t, e in zip(r["tensiones"], (6 * _r2, 6 * _r2, -2 * _r2))] == [0, 0, 0]
                        and not r["todas_positivas"] and not r["sostiene"] and r["equilibrio"]
                        and len(r["avisos"]) == 1 and r["avisos"][0].startswith("el cable 3 tendría que empujar"))(
         equilibrio_particula([(1, 0, 1), (-1, 0, 1), (0, 1, 1)], (0, 2, -10)))),
    ("anclajes coplanares (1, 0, 1), (-1, 0, 1), (2, 0, 1) → ValueError «los cables están en un mismo plano»",
     lambda: _mensaje(lambda: equilibrio_particula("1, 0, 1; -1, 0, 1; 2, 0, 1", "0, 0, -10")).startswith(
         "los cables están en un mismo plano: con esta geometría no hay solución única")),
    ("momento((0, 0, -200), (0.3, 0.4, 0)) = (-80, 60, 0), con M·r = 0 y M·F = 0",
     lambda: (lambda M, r, F: M == vector("-80, 60, 0") and _punto(M, r) == 0 and _punto(M, F) == 0)(
         momento((0, 0, -200), (0.3, 0.4, 0)), vector("0.3, 0.4, 0"), vector("0, 0, -200"))),
    ("momento respecto a O = (1, 0, 0) de F = (0, 10, 0) aplicada en P = (3, 0, 0) → (0, 0, 20) (y no el del origen, "
     "(0, 0, 30))",
     lambda: momento((0, 10, 0), (3, 0, 0), (1, 0, 0)) == vector("0, 0, 20")
             and momento((0, 10, 0), (3, 0, 0)) == vector("0, 0, 30")),
    ("trabajo((40, 30, 0), (5, 0, 0)) = 200 y con d = (3, 4, 0) = 240; con una fuerza perpendicular a d, 0",
     lambda: trabajo((40, 30, 0), (5, 0, 0)) == 200 and trabajo((40, 30, 0), (3, 4, 0)) == 240
             and trabajo((40, 30, 0), (-3, 4, 7)) == 0),
    ("vector_unitario((1, 1, 1), (1, 1, 1)) → ValueError «P y Q son el mismo punto»",
     lambda: _mensaje(lambda: vector_unitario((1, 1, 1), (1, 1, 1))).startswith("P y Q son el mismo punto")),
    ("prueba de humo: la calculadora corre con sus valores iniciales y con un cable que empujaría, sin error y sin "
     "«Revisa la entrada»",
     lambda: _humo(calculadora_equilibrio, _CONTROLES_1)
             and _humo(calculadora_equilibrio, _CONTROLES_1, anclajes_txt="1, 0, 1\n-1, 0, 1\n0, 1, 1",
                       carga_txt="0, 2, -10", unidad="kN", n_cifras=6)),
]
corre_pruebas(PRUEBAS_1, 1)

# %% [markdown]
# **Contrasta (sección 1).** Calcula a mano el momento respecto al origen de $\mathbf F=(0,50,0)$ N aplicada en $(0.2,\,0,\,0)$ m: primero $\mathbf r$, después $\mathbf r\times\mathbf F$ componente por componente. Escribe `mi_momento = (Mx, My, Mz)` en N·m.

# %%
mi_momento = None     # tres números (Mx, My, Mz) en N·m, por ejemplo (1, -2, 0.5)

if mi_momento is None:     # la referencia se muestra después de tu cálculo
    print("Falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda.")
else:
    try:
        mio = vector(mi_momento, "mi_momento")
    except ValueError as err:
        print(f"Revisa tu respuesta: {str(err).rstrip('.')}.")
    else:
        if mio.rows != 3:          # respuesta imposible: la referencia todavía no se muestra
            print(f"El momento es un vector de R³, con tres componentes (escribiste {mio.rows}): revisa tu respuesta y "
                  "vuelve a ejecutar la celda.")
        else:
            F_c, P_c = vector("0, 50, 0"), vector("0.2, 0, 0")          # los datos del enunciado
            ref = momento(F_c, P_c)
            print(f"La calculadora da: M = r × F = {_tupla_u(ref, 'N·m')}.")
            if all(abs(float(a - b)) <= 0.01 for a, b in zip(mio, ref)):
                print("Coinciden (con tolerancia de 0.01 N·m por componente).")
            elif all(abs(float(a + b)) <= 0.01 for a, b in zip(mio, ref)):
                print("NO coinciden: tu vector es el opuesto. Probablemente calculaste F × r; el momento es r × F, con r "
                      "del origen al punto de aplicación (el orden importa: F × r = -M).")
            else:
                print("NO coinciden: escribe r = (0.2, 0, 0) y F = (0, 50, 0) y aplica r × F = (ry·Fz - rz·Fy, "
                      "rz·Fx - rx·Fz, rx·Fy - ry·Fx).")
            print("\nEl procedimiento completo:")
            _informe_momento(F_c, P_c, (0, 0, 0), 4)
