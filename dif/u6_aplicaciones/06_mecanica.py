# ID: DIF-U6-NB06
# Notebook: dif/u6_aplicaciones.ipynb · sección 6.6 mecánica
# Repositorio: dif/u6_aplicaciones/06_mecanica.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6. Mecánica: velocidad, aceleración y dimensiones óptimas
#
# $v=s'$ y $a=v'=s''$. En reposo donde $v=0$; la distancia recorrida suma los tramos entre cambios de sentido; el desplazamiento es $s(t_1)-s(t_0)$.
#
# Dos calculadoras: la primera analiza un movimiento $s(t)$ (escribe $t$ como `x`); la segunda dimensiona un tanque cilíndrico cerrado de volumen $V$ con precios distintos para tapas y pared y un radio máximo.

# %%
def movimiento(s_txt, t0, t1):
    """Diccionario con v, a, instantes de reposo, distancia recorrida y desplazamiento en [t0, t1]."""
    s = parsear(s_txt)
    t0, t1 = sp.nsimplify(t0), sp.nsimplify(t1)
    if not t0 < t1:
        raise ValueError("el intervalo debe cumplir t0 < t1.")
    v, a = sp.diff(s, x), sp.diff(s, x, 2)
    reposo = sorted(sp.solveset(v, x, domain=sp.Interval.open(t0, t1)), key=float)
    marcas = [t0, *reposo, t1]
    distancia = sum(abs(s.subs(x, q) - s.subs(x, p)) for p, q in zip(marcas[:-1], marcas[1:]))
    return {"v": v, "a": a, "reposo": reposo, "distancia": sp.simplify(distancia),
            "desplazamiento": sp.simplify(s.subs(x, t1) - s.subs(x, t0))}


def tanque_optimo(V, precio_tapas, precio_pared, r_max=None):
    """(r, h, costo) del tanque cilíndrico cerrado de costo mínimo; respeta r <= r_max si se da."""
    if min(V, precio_tapas, precio_pared) <= 0:
        raise ValueError("volumen y precios deben ser positivos.")
    r = sp.symbols("r", positive=True)
    C = precio_tapas * 2 * sp.pi * r**2 + precio_pared * 2 * V / r
    ro = sp.solve(sp.diff(C, r), r)[0]
    if r_max is not None and ro > r_max:
        ro = sp.nsimplify(r_max)          # C decrece antes del óptimo: el mejor radio permitido es el tope
    return float(ro), float(V / (sp.pi * ro**2)), float(C.subs(r, ro))


def calculadora_movimiento(s_txt, t0, t1, n_cifras):
    try:
        r = movimiento(s_txt, t0, t1)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print("v(t) =", r["v"], "   a(t) =", r["a"])
    print("reposo en t =", [cifras(q, n_cifras) for q in r["reposo"]] or "ningún instante del intervalo")
    print(f"distancia recorrida = {cifras(r['distancia'], n_cifras)}   desplazamiento = {cifras(r['desplazamiento'], n_cifras)}")


widgets.interact(calculadora_movimiento,
    s_txt=widgets.Text(value="20*x - 2.5*x**2", description="s(t) ="),
    t0=widgets.FloatText(value=0.0, description="t0"), t1=widgets.FloatText(value=6.0, description="t1"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));


def calculadora_tanque(V, precio_tapas, precio_pared, r_max, n_cifras):
    try:
        r, h, c = tanque_optimo(V, precio_tapas, precio_pared, r_max if r_max > 0 else None)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"r = {cifras(r, n_cifras)} m, h = {cifras(h, n_cifras)} m, h/r = {cifras(h / r, n_cifras)},"
          f" costo = {cifras(c, n_cifras)} pesos  ({n_cifras} cifras significativas)")


widgets.interact(calculadora_tanque,
    V=widgets.FloatText(value=0.5, description="V (m³)"),
    precio_tapas=widgets.FloatText(value=675, description="tapas $/m²"),
    precio_pared=widgets.FloatText(value=450, description="pared $/m²"),
    r_max=widgets.FloatText(value=0, description="r máx (0 = sin tope)"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 6 (resultado conocido)
_q = movimiento("0.4*(10*(x/2)**3 - 15*(x/2)**4 + 6*(x/2)**5)", 0, 2)
PRUEBAS_6 = [
    ("frenado 20t - 2.5t²: reposo en 4 s, 40 m", lambda: movimiento("20*x - 2.5*x**2", 0, 6)["reposo"] == [4]
                                                         and movimiento("20*x - 2.5*x**2", 0, 4)["desplazamiento"] == 40),
    ("t³ - 9t² + 24t en [0, 5]: reposo 2 y 4; distancia 28; desplazamiento 20",
     lambda: movimiento("x**3 - 9*x**2 + 24*x", 0, 5)["reposo"] == [2, 4]
             and movimiento("x**3 - 9*x**2 + 24*x", 0, 5)["distancia"] == 28
             and movimiento("x**3 - 9*x**2 + 24*x", 0, 5)["desplazamiento"] == 20),
    ("quíntico: v máx 0.375 en t = 1", lambda: _q["v"].subs(x, 1) == sp.Rational(3, 8)),
    ("quíntico: v y a nulas en 0 y 2", lambda: [_q["v"].subs(x, k) for k in (0, 2)] == [0, 0] and [_q["a"].subs(x, k) for k in (0, 2)] == [0, 0]),
    ("tanque 0.5 m³, 675/450: r ≈ 0.3758, h/r = 3, ≈ 1796.4",
     lambda: cerca(tanque_optimo(0.5, 675, 450)[0], 0.375751, 1e-5) and cerca(tanque_optimo(0.5, 675, 450)[1] / tanque_optimo(0.5, 675, 450)[0], 3, 1e-9)
             and round(tanque_optimo(0.5, 675, 450)[2], 1) == 1796.4),
    ("tanque con r ≤ 0.35: ≈ 1805.26", lambda: round(tanque_optimo(0.5, 675, 450, 0.35)[2], 2) == 1805.26),
    ("precios iguales: h = 2r", lambda: cerca(tanque_optimo(1.0, 100, 100)[1] / tanque_optimo(1.0, 100, 100)[0], 2, 1e-9)),
]
for nombre, prueba in PRUEBAS_6:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 6).** Para $s(t)=2t^3-9t^2+12t$ (m, s), encuentra a mano los instantes de reposo en $[0,3]$ y escríbelos como lista.

# %%
mis_instantes = None     # escribe una lista, por ejemplo: [0.5, 2.5]

ref = movimiento("2*x**3 - 9*x**2 + 12*x", 0, 3)["reposo"]
if mis_instantes is None:                     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    print("La calculadora da:", ref)
    print("coinciden" if sorted(mis_instantes) == ref else "NO coinciden: v = 6t² - 18t + 12 = 6(t - 1)(t - 2)")
