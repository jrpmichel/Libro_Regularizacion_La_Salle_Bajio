# ID: ALG-U1-NB03
# Notebook: alg/u1_vectores.ipynb · sección 1.3 producto punto y proyección
# Repositorio: alg/u1_vectores/03_punto.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Producto punto y proyección
#
# El producto punto de dos vectores con la misma cantidad de componentes es un **número**:
#
# $$\mathbf a\cdot\mathbf b=a_1b_1+a_2b_2+\dots+a_nb_n=\|\mathbf a\|\,\|\mathbf b\|\cos\theta,\qquad 0\le\theta\le\pi .$$
#
# * **Ángulo:** $\theta=\arccos\dfrac{\mathbf a\cdot\mathbf b}{\|\mathbf a\|\,\|\mathbf b\|}$, con $\mathbf a,\mathbf b\neq\mathbf 0$. Por redondeo, el cociente de dos vectores paralelos puede salir $1.0000000000000002$, y `arccos` de un número mayor que 1 da `nan`. Por eso la calculadora lo **recorta** al intervalo $[-1,1]$ (`np.clip`) antes de aplicar `arccos`.
# * **Signo:** $\mathbf a\cdot\mathbf b>0$, ángulo agudo; $\mathbf a\cdot\mathbf b=0$, recto (los vectores son **ortogonales**); $\mathbf a\cdot\mathbf b<0$, obtuso.
# * **Componente y proyección** de $\mathbf a$ sobre $\mathbf b\neq\mathbf 0$, y la parte perpendicular:
#
# $$\mathrm{comp}_{\mathbf b}\,\mathbf a=\frac{\mathbf a\cdot\mathbf b}{\|\mathbf b\|},\qquad \mathrm{proy}_{\mathbf b}\,\mathbf a=\frac{\mathbf a\cdot\mathbf b}{\mathbf b\cdot\mathbf b}\,\mathbf b,\qquad \mathbf a_\perp=\mathbf a-\mathrm{proy}_{\mathbf b}\,\mathbf a .$$
#
# La componente es negativa si el ángulo es obtuso. La comprobación es $\mathbf a_\perp\cdot\mathbf b=0$, y $\mathbf a=\mathrm{proy}_{\mathbf b}\,\mathbf a+\mathbf a_\perp$.
#
# **Similitud coseno.** En ciencia de datos, $\cos\theta$ entre dos vectores de $n$ componentes (por ejemplo, cuántas veces aparece cada palabra en dos documentos) mide qué tan parecidas son sus direcciones sin importar su tamaño: 1 si apuntan igual, 0 si son ortogonales. Por eso la calculadora acepta vectores de 2, 3 o más componentes.

# %%
def punto(a, b):
    """a · b = a1 b1 + a2 b2 + ... (misma cantidad de componentes)."""
    a, b = _par(a, b)
    return _punto(a, b)


def _coseno(a, b):
    """cos θ entre a y b ya leídos, recortado a [-1, 1]. ValueError si alguno es el vector cero."""
    if not np.any(a) or not np.any(b):
        raise ValueError("el vector cero no tiene dirección: no forma ángulo con ningún vector.")
    ea, eb = _escalado(a), _escalado(b)       # escalar por potencias de 2 evita desbordes y no cambia θ
    c = float(np.clip(_punto(ea, eb) / (_norma(ea) * _norma(eb)), -1.0, 1.0))
    # arccos no distingue ángulos menores que ~1e-8 rad: un coseno a 1e-15 de ±1 es de vectores paralelos
    return math.copysign(1.0, c) if abs(c) > 1 - 1e-15 else c


def angulo(a, b, unidad="grados"):
    """Ángulo entre a y b, de 0 a 180° (o de 0 a π rad), con arccos del coseno recortado."""
    a, b = _par(a, b)
    return en_unidad(math.acos(_coseno(a, b)), unidad)


def clasifica(a, b, tol=1e-12):
    """'agudo', 'recto' u 'obtuso' según el signo de a · b (con |cos θ| <= tol se toma como recto)."""
    a, b = _par(a, b)
    c = _coseno(a, b)
    return "recto" if abs(c) <= tol else ("agudo" if c > 0 else "obtuso")


def proyeccion(a, b):
    """(comp_b a, proy_b a, a⊥). b no puede ser el vector cero."""
    a, b = _par(a, b)
    if not np.any(b):
        raise ValueError("no se proyecta sobre el vector cero: b no tiene dirección.")
    ub = b / _norma(b)                        # (a·b/‖b‖) b/‖b‖ = (a·b / b·b) b, sin desbordes
    comp = _punto(a, ub)
    par = comp * ub
    return comp, par, _limpia(a - par, np.abs(a) + np.abs(par))


def calculadora_punto(a_txt, b_txt, unidad, n_cifras):
    try:
        a, b = _par(a_txt, b_txt)
        p = punto(a, b)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"a = {fmt(a, n_cifras)}    b = {fmt(b, n_cifras)}")
    print(f"a · b = {cifras(p, n_cifras)}")
    try:
        t, clase = angulo(a, b, unidad), clasifica(a, b)
        print(f"cos θ = {cifras(_coseno(a, b), n_cifras)}   (la similitud coseno)")
        signo = {"agudo": "a · b > 0", "obtuso": "a · b < 0",
                 "recto": "a · b = 0" if p == 0 else "a · b ≈ 0"}[clase]
        print(f"θ = {cifras(t, n_cifras)}{_sufijo(unidad)}   →   ángulo {clase} ({signo})")
    except ValueError as err:
        print("ángulo:", err)
    try:
        comp, par, perp = proyeccion(a, b)
    except ValueError as err:
        print("proyección:", err); return
    print(f"comp_b a = {cifras(comp, n_cifras)}")
    print(f"proy_b a = {fmt(par, n_cifras)}   (la parte de a paralela a b)")
    print(f"a⊥ = a − proy_b a = {fmt(perp, n_cifras)}   (la parte de a perpendicular a b)")
    r = punto(perp, b)
    if abs(r) <= 1e-9 * _norma(a) * _norma(b):
        print(f"comprobación: a⊥ · b = {cifras(r, n_cifras) if r else '0'}{' ≈ 0' if r else ''}, "
              "así que a⊥ es perpendicular a b")
    else:
        print(f"comprobación: a⊥ · b = {cifras(r, n_cifras)} ≠ 0: revisa la entrada")


widgets.interact(calculadora_punto,
    a_txt=widgets.Text(value="4, 1, 1", description="a =", continuous_update=False),
    b_txt=widgets.Text(value="2, 2, 1", description="b =", continuous_update=False),
    unidad=widgets.Dropdown(options=["grados", "radianes"], value="grados", description="unidad"),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras sig.", continuous_update=False));

# %%
# Casos de prueba de la sección 3 (resultado conocido)
PRUEBAS_3 = [
    ("(4, 1, 1) · (2, 2, 1) = 11, ángulo 30.20° = 0.5272 rad",
     lambda: cerca(punto((4, 1, 1), (2, 2, 1)), 11) and round(angulo((4, 1, 1), (2, 2, 1), "grados"), 2) == 30.20
             and round(angulo((4, 1, 1), (2, 2, 1), "radianes"), 4) == 0.5272),
    ("clasificación: (1, 2, -1) y (3, 0, 3) recto; (4, 1, 1) y (2, 2, 1) agudo; (1, 1, 0) y (-1, 0, 1) obtuso",
     lambda: clasifica((1, 2, -1), (3, 0, 3)) == "recto" and clasifica((4, 1, 1), (2, 2, 1)) == "agudo"
             and clasifica((1, 1, 0), (-1, 0, 1)) == "obtuso"),
    ("proyección de (3, 1, 2) sobre (1, 2, 2): componente 3, paralela (1, 2, 2), perpendicular (2, -1, 0)",
     lambda: cerca(proyeccion((3, 1, 2), (1, 2, 2))[0], 3)
             and np.allclose(proyeccion((3, 1, 2), (1, 2, 2))[1], (1, 2, 2))
             and np.allclose(proyeccion((3, 1, 2), (1, 2, 2))[2], (2, -1, 0))),
    ("ángulo entre (1, 0) y (0, 1) = π/2 rad", lambda: cerca(angulo((1, 0), (0, 1), "radianes"), math.pi / 2)),
    ("ángulo entre (1, 1, 1) y (2, 2, 2) = 0 (sin nan)", lambda: angulo((1, 1, 1), (2, 2, 2)) == 0),
    ("el ángulo con el vector cero se rechaza", lambda: _rechaza(lambda: angulo((0, 0, 0), (1, 2, 3)))),
]
for nombre, prueba in PRUEBAS_3:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 3).** Antes de calcular, **predice** con el signo de $(2,1,-2)\cdot(3,0,4)$ si el ángulo entre los dos vectores es agudo, recto u obtuso. Después calcula a mano el ángulo y escríbelo en grados o en radianes (indica cuál). Compara con la calculadora en las dos unidades.

# %%
mi_prediccion = None   # "agudo", "recto" u "obtuso", según el signo de a · b
mi_angulo = None       # tu ángulo, por ejemplo: 45.0 o "pi/4"
mi_unidad = "grados"   # "grados" o "radianes": la unidad en que escribiste tu ángulo

if mi_prediccion is None or mi_angulo is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    a_c, b_c = (2, 1, -2), (3, 0, 4)
    print(f"La calculadora da: a · b = {cifras(punto(a_c, b_c), 4)}, ángulo {clasifica(a_c, b_c)}, "
          f"θ = {cifras(angulo(a_c, b_c, 'grados'), 5)}° = {cifras(angulo(a_c, b_c, 'radianes'), 5)} rad")
    pred = str(mi_prediccion).strip(" .").lower()
    print("predicción: correcta" if pred == clasifica(a_c, b_c) else
          "predicción: NO coincide; el signo de a · b decide: positivo agudo, cero recto, negativo obtuso")
    unidad_c = str(mi_unidad).strip().lower()
    otra = {"grados": "radianes", "radianes": "grados"}
    try:
        texto = mi_angulo.replace("°", "").replace("º", "") if isinstance(mi_angulo, str) else mi_angulo
        valor = numero(texto, "tu ángulo")
        ok = cerca(math.degrees(de_unidad(valor, unidad_c)), angulo(a_c, b_c, "grados"), 0.05)   # tolerancia: 0.05°
        if not ok and cerca(math.degrees(de_unidad(valor, otra[unidad_c])), angulo(a_c, b_c, "grados"), 0.05):
            print(f"Tu número está bien, pero en {otra[unidad_c]}: cambia mi_unidad a \"{otra[unidad_c]}\".")
    except ValueError as err:
        print("Revisa tu valor:", err); ok = None
    if ok is not None:
        print("ángulo: coinciden" if ok else
              "ángulo: NO coinciden; θ = arccos(a · b / (‖a‖ ‖b‖)), con la calculadora en la unidad que escribiste")
