# ID: DIF-U4-NB05
# Notebook: dif/u4_derivada.ipynb · sección 4.5 la derivada como razón de cambio
# Repositorio: dif/u4_derivada/05_razon_cambio.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 5. La derivada como razón de cambio
#
# La pendiente de la secante entre $a$ y $a+h$ es la **razón de cambio promedio**; la derivada $f'(a)$ es la **razón de cambio instantánea**. Las dos tienen las mismas unidades: las de $y$ entre las de $x$ (m/s, °C/min, pesos/pieza, A/s).
#
# La calculadora recibe $f$ (con $x$ como variable, aunque represente el tiempo), el instante $a$, el intervalo $h$ y las unidades. Da las dos razones de cambio con sus unidades y cuánto se separan.

# %%
def razones(f, a, hh):
    """(promedio en [a, a + h], instantánea en a), exactas."""
    a, hh = sp.nsimplify(a), sp.nsimplify(hh)
    prom = pendiente_secante(f, a, a + hh)
    izq, der, inst = derivada_en(f, a)
    if inst is None:
        raise ValueError(f"la razón instantánea no existe en x = {a} (izquierda: {lado_texto(izq)}; derecha: {lado_texto(der)}).")
    return prom, inst


def calculadora_razon(f_txt, a, h_val, u_y, u_x, n_cifras):
    try:
        f = parsear(f_txt)
        prom, inst = razones(f, a, h_val)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    u = f"{u_y}/{u_x}" if u_y and u_x else "(sin unidades)"
    print(f"promedio en [{a}, {a + h_val}]: {cifras(prom, n_cifras)} {u}")
    print(f"instantánea en {a}: {cifras(inst, n_cifras)} {u}")
    dif = sp.N(prom - inst)
    rel = f" ({cifras(abs(dif / inst) * 100, 3)} % de la instantánea)" if inst != 0 else ""
    print(f"diferencia: {cifras(dif, 3)} {u}{rel}   ({n_cifras} cifras significativas)")


widgets.interact(calculadora_razon,
    f_txt=widgets.Text(value="0.6*x**2 - 0.1*x**3", description="f(x) ="),
    a=widgets.FloatText(value=2.0, description="a"), h_val=widgets.FloatText(value=1.0, description="h"),
    u_y=widgets.Text(value="m", description="unid. y"), u_x=widgets.Text(value="s", description="unid. x"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 5 (resultado conocido)
PRUEBAS_5 = [
    ("caída libre 4.9t² en 2: promedio con h = 1 es 24.5 m/s, instantánea 19.6 m/s",
     lambda: razones(parsear("4.9*x**2"), 2, 1) == (sp.Rational(49, 2), sp.Rational(98, 5))),
    ("eje del robot en 2: promedio con h = 1 es 1.1 m/s, instantánea 1.2 m/s",
     lambda: razones(parsear("0.6*x**2 - 0.1*x**3"), 2, 1) == (sp.Rational(11, 10), sp.Rational(6, 5))),
    ("costo en q = 100: promedio con h = 20 es 26 pesos/pieza, instantánea 25",
     lambda: razones(parsear("2000 + 15*x + 0.05*x**2"), 100, 20) == (26, 25)),
    ("con h negativo se obtiene el promedio por la izquierda: 4.9t² en 2 con h = -1 da 14.7",
     lambda: razones(parsear("4.9*x**2"), 2, -1)[0] == sp.Rational(147, 10)),
    ("|x| en 0 no tiene razón instantánea", lambda: _rechaza(lambda: razones(sp.Abs(x), 0, 1))),
]
for nombre, prueba in PRUEBAS_5:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 5).** Para $V(t)=3t^2+2$, calcula a mano la razón de cambio promedio entre $t=1$ y $t=2$ y la instantánea en $t=1$. Compara con la calculadora ($a=1$, $h=1$).

# %%
mi_promedio, mi_instantanea = None, None     # escribe dos números, por ejemplo: 4, 2

ref = razones(3*x**2 + 2, 1, 1)
print("La calculadora da: promedio", ref[0], "| instantánea", ref[1])
print("falta tu cálculo" if mi_promedio is None else
      ("coincide" if cerca(mi_promedio, float(ref[0])) and cerca(mi_instantanea, float(ref[1]))
       else "NO coincide: el promedio es (V(2) - V(1))/1 y la instantánea es V'(1) = 6·1"))
