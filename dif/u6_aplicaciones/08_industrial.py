# ID: DIF-U6-NB08
# Notebook: dif/u6_aplicaciones.ipynb · sección 6.8 industrial
# Repositorio: dif/u6_aplicaciones/08_industrial.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 8. Industrial: utilidad máxima y lote económico
#
# **Utilidad.** $U=R-C$ con $R=q\,p(q)$. En un máximo interior, $R'(q)=C'(q)$ (ingreso marginal igual a costo marginal) y $U''<0$. Si hay capacidad máxima, se compara con el extremo.
#
# **Lote económico.** $CT(Q)=\dfrac{DS}{Q}+\dfrac{HQ}{2}$, mínimo en $Q^*=\sqrt{2DS/H}$. Supuestos: demanda constante, reposición instantánea, sin faltantes ni descuentos.
#
# Escribe la cantidad $q$ como `x`.

# %%
def utilidad_maxima(p_txt, C_txt, q_max):
    """(q de utilidad máxima en [0, q_max], U, q de ingreso máximo, U en ese q)."""
    p, C = parsear(p_txt), parsear(C_txt)
    R, U = x * p, x * p - C
    qU = [c for c in sp.solveset(sp.diff(U, x), x, domain=sp.Interval.open(0, q_max))]
    cand = [sp.Integer(0), sp.nsimplify(q_max), *qU]
    mejor = max(cand, key=lambda c: float(U.subs(x, c)))
    qR = [c for c in sp.solveset(sp.diff(R, x), x, domain=sp.Interval.open(0, sp.oo))]
    qR = qR[0] if qR else None
    return mejor, U.subs(x, mejor), qR, (U.subs(x, qR) if qR is not None else None)


def lote_economico(D, S, H, Q_max=None):
    """(Q*, costo anual en Q*, Q usado, costo con Q usado)."""
    if min(D, S, H) <= 0:
        raise ValueError("D, S y H deben ser positivos.")
    Qs = math.sqrt(2 * D * S / H)
    CT = lambda Q: D * S / Q + H * Q / 2
    Qu = min(Qs, Q_max) if Q_max else Qs        # CT decrece antes de Q*: con tope, el mejor es el tope
    return Qs, CT(Qs), Qu, CT(Qu)


def calculadora_utilidad(p_txt, C_txt, q_max, n_cifras):
    try:
        q, U, qR, UR = utilidad_maxima(p_txt, C_txt, q_max)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"utilidad máxima en [0, {q_max}]: q = {cifras(q, n_cifras)}, U = {cifras(U, n_cifras)}  ({n_cifras} cifras significativas)")
    if qR is not None:
        print(f"(el ingreso es máximo en q = {cifras(qR, n_cifras)}, donde U = {cifras(UR, n_cifras)})")


widgets.interact(calculadora_utilidad,
    p_txt=widgets.Text(value="500 - 0.5*x", description="p(q) ="),
    C_txt=widgets.Text(value="20000 + 140*x", description="C(q) ="),
    q_max=widgets.FloatText(value=1000, description="capacidad"),
    n_cifras=widgets.IntSlider(value=6, min=1, max=8, description="cifras sig."));


def calculadora_lote(D, S, H, Q_max, n_cifras):
    try:
        Qs, CTs, Qu, CTu = lote_economico(D, S, H, Q_max if Q_max > 0 else None)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"Q* = {cifras(Qs, n_cifras)}, costo anual {cifras(CTs, n_cifras)}  ({n_cifras} cifras significativas)")
    if Qu != Qs:
        print(f"con el tope: Q = {cifras(Qu, n_cifras)}, costo anual {cifras(CTu, n_cifras)}"
              f" ({cifras(100 * (CTu / CTs - 1), 3)} % más)")


widgets.interact(calculadora_lote,
    D=widgets.FloatText(value=12000, description="D (unid/año)"), S=widgets.FloatText(value=800, description="S ($/pedido)"),
    H=widgets.FloatText(value=30, description="H ($/unid·año)"), Q_max=widgets.FloatText(value=0, description="tope (0 = no)"),
    n_cifras=widgets.IntSlider(value=6, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 8 (resultado conocido)
PRUEBAS_8 = [
    ("utilidad: q = 360, U = 44800", lambda: utilidad_maxima("500 - 0.5*x", "20000 + 140*x", 1000)[:2] == (360, 44800)),
    ("ingreso máximo en 500 con U = 35000", lambda: utilidad_maxima("500 - 0.5*x", "20000 + 140*x", 1000)[2:] == (500, 35000)),
    ("con capacidad 300: q = 300, U = 43000", lambda: utilidad_maxima("500 - 0.5*x", "20000 + 140*x", 300)[:2] == (300, 43000)),
    ("EOQ 12000/800/30: Q* = 800, 24000", lambda: lote_economico(12000, 800, 30)[:2] == (800.0, 24000.0)),
    ("EOQ con tope 600: 25000", lambda: lote_economico(12000, 800, 30, 600)[2:] == (600, 25000.0)),
    ("EOQ proveedor B 12000/1100/30: Q* ≈ 938.1, ≈ 28142.5", lambda: round(lote_economico(12000, 1100, 30)[0], 1) == 938.1
                                                                  and round(lote_economico(12000, 1100, 30)[1], 1) == 28142.5),
    ("datos no positivos se rechazan", lambda: _rechaza(lambda: lote_economico(0, 800, 30))),
]
for nombre, prueba in PRUEBAS_8:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 8).** Con $D=5000$, $S=400$ y $H=25$, calcula a mano $Q^*$ y escribe el número.

# %%
mi_Q = None              # escribe un número

ref = lote_economico(5000, 400, 25)[0]
if mi_Q is None:                     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    print("La calculadora da Q* =", cifras(ref, 6))
    print("coincide" if abs(mi_Q - ref) < 0.5 else "NO coincide: Q* = raíz de 2·D·S/H")
