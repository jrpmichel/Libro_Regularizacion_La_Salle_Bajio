# ID: DIF-U6-NB07
# Notebook: dif/u6_aplicaciones.ipynb · sección 6.7 eléctrica y señales
# Repositorio: dif/u6_aplicaciones/07_electrica.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 7. Eléctrica y señales: $i_C=C\,dv/dt$, $v_L=L\,di/dt$; máximos de una señal
#
# Dos calculadoras: la de capacitor da $i_C$ a partir de $C$ y $v(t)$; la de inductor, $v_L$ a partir de $L$ e $i(t)$. Las dos buscan el máximo y el mínimo del resultado en un intervalo. Convención pasiva: $i$ entra por la terminal $+$.
#
# Unidades: $C$ en faradios, $L$ en henrios, $v$ en volts, $i$ en amperes, $t$ en segundos (escribe $t$ como `x`). Los extremos de una señal en un intervalo están donde su derivada se anula **o en los extremos del intervalo**.

# %%
def corriente_capacitor(C, v_txt):
    if C <= 0:
        raise ValueError("la capacitancia debe ser positiva (en faradios: 10 µF = 10e-6).")
    return sp.nsimplify(C) * sp.diff(parsear(v_txt), x)


def tension_inductor(L, i_txt):
    if L <= 0:
        raise ValueError("la inductancia debe ser positiva (en henrios: 20 mH = 0.02).")
    return sp.nsimplify(L) * sp.diff(parsear(i_txt), x)


def extremos_senal(senal, t0, t1, n=4001):
    """(t del máximo, valor, t del mínimo, valor) en [t0, t1]: puntos críticos refinados y extremos del intervalo."""
    f = sp.lambdify(x, senal, "numpy")
    ts = np.linspace(float(t0), float(t1), n)
    with np.errstate(all="ignore"):
        vs = np.broadcast_to(f(ts), ts.shape).astype(float)
    cand = [float(t0), float(t1)]
    dfun = sp.lambdify(x, sp.diff(senal, x), "numpy")
    signos = np.sign(np.broadcast_to(dfun(ts), ts.shape))
    for k in np.where(signos[:-1] * signos[1:] < 0)[0]:    # cambios de signo de la derivada
        a, b = ts[k], ts[k + 1]
        for _ in range(60):                                # bisección sobre la derivada
            m = (a + b) / 2
            a, b = (a, m) if np.sign(dfun(a)) * np.sign(dfun(m)) <= 0 else (m, b)
        cand.append(float((a + b) / 2))
    vals = [float(senal.subs(x, c)) for c in cand]
    kM, km = int(np.argmax(vals)), int(np.argmin(vals))
    return cand[kM], vals[kM], cand[km], vals[km]


def _reporte(nombre, unidad, senal, t0_ms, t1_ms, n):
    tM, sM, tm, sm = extremos_senal(senal, t0_ms / 1000, t1_ms / 1000)
    print(f"{nombre}(t) =", sp.N(senal, n), unidad)
    print(f"en [{t0_ms}, {t1_ms}] ms: máximo {cifras(sM, n)} {unidad} en t = {cifras(tM * 1000, n)} ms"
          f" | mínimo {cifras(sm, n)} {unidad} en t = {cifras(tm * 1000, n)} ms  ({n} cifras significativas)")


def calculadora_capacitor(C_uF, v_txt, t0_ms, t1_ms, n_cifras):
    try:
        _reporte("i_C", "A", corriente_capacitor(C_uF * 1e-6, v_txt), t0_ms, t1_ms, n_cifras)
    except (ValueError, TypeError) as err:
        print("Revisa la entrada:", err)


widgets.interact(calculadora_capacitor,
    C_uF=widgets.FloatText(value=2.0, description="C (µF)"),
    v_txt=widgets.Text(value="50*(exp(-500*x) - exp(-2000*x))", description="v(t) ="),
    t0_ms=widgets.FloatText(value=0.0, description="t0 (ms)"), t1_ms=widgets.FloatText(value=6.0, description="t1 (ms)"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));


def calculadora_inductor(L_mH, i_txt, t0_ms, t1_ms, n_cifras):
    try:
        _reporte("v_L", "V", tension_inductor(L_mH * 1e-3, i_txt), t0_ms, t1_ms, n_cifras)
    except (ValueError, TypeError) as err:
        print("Revisa la entrada:", err)


widgets.interact(calculadora_inductor,
    L_mH=widgets.FloatText(value=20.0, description="L (mH)"),
    i_txt=widgets.Text(value="3*sin(377*x)", description="i(t) ="),
    t0_ms=widgets.FloatText(value=0.0, description="t0 (ms)"), t1_ms=widgets.FloatText(value=20.0, description="t1 (ms)"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 7 (resultado conocido)
_pulso = parsear("50*(exp(-500*x) - exp(-2000*x))")
PRUEBAS_7 = [
    ("10 µF con 10 sen 377t: pico 0.0377 A", lambda: cerca(corriente_capacitor(10e-6, "10*sin(377*x)").subs(x, 0), 0.0377, 1e-12)),
    ("20 mH con 3 sen 377t: pico 22.62 V", lambda: cerca(tension_inductor(0.02, "3*sin(377*x)").subs(x, 0), 22.62, 1e-9)),
    ("pulso: v máx ≈ 23.62 V en ≈ 0.924 ms", lambda: round(extremos_senal(_pulso, 0, 0.006)[1], 2) == 23.62
                                                    and round(extremos_senal(_pulso, 0, 0.006)[0] * 1000, 3) == 0.924),
    ("pulso: i máx 0.15 A en t = 0 (extremo del intervalo)",
     lambda: extremos_senal(corriente_capacitor(2e-6, "50*(exp(-500*x) - exp(-2000*x))"), 0, 0.006)[:2] == (0.0, 0.15)),
    ("pulso: i mín ≈ -0.0149 A en ≈ 1.85 ms",
     lambda: round(extremos_senal(corriente_capacitor(2e-6, "50*(exp(-500*x) - exp(-2000*x))"), 0, 0.006)[3], 4) == -0.0149),
    ("capacitancia no positiva se rechaza", lambda: _rechaza(lambda: corriente_capacitor(0, "x"))),
]
for nombre, prueba in PRUEBAS_7:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 7).** Calcula a mano la corriente pico de un capacitor de 47 µF con $v=5\,\mathrm{sen}\,1000t$ y escribe el número en amperes.

# %%
mi_pico = None           # escribe un número

ref = float(corriente_capacitor(47e-6, "5*sin(1000*x)").subs(x, 0))
if mi_pico is None:                     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    print("La calculadora da:", cifras(ref, 4), "A")
    print("coincide" if abs(mi_pico - ref) < 5e-4 else "NO coincide: i = C·5·1000·cos(1000t) y el pico es el coeficiente")
