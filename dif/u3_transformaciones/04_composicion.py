# ID: DIF-U3-NB04
# Notebook: dif/u3_transformaciones.ipynb · sección 3.4 composición de transformaciones
# Repositorio: dif/u3_transformaciones/04_composicion.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 4. Composición de transformaciones
#
# Si la fórmula trae $f(bx+c)$, primero se **factoriza** $b$: $f(bx+c)=f\big(b(x-h)\big)$ con $h=-c/b$. El desplazamiento es $h$, no $c$.
#
# La calculadora recibe $b$ y $c$, da $h$ y compara la función correcta con la que resulta de desplazar $c$ en lugar de $h$. También muestra los dos órdenes posibles de los pasos horizontales.

# %%
def factorizar_interior(b, c):
    """h tal que b x + c = b (x - h)."""
    b, c = sp.nsimplify(b), sp.nsimplify(c)
    if b == 0:
        raise ValueError("b = 0: el interior no depende de x.")
    return -c / b


def calculadora_orden(base, b, c):
    try:
        f = BASES[base][0]
        h = factorizar_interior(b, c)
    except (KeyError, ValueError) as err:
        print("Revisa la entrada:", err); return
    b_, c_ = sp.nsimplify(b), sp.nsimplify(c)
    lado = lambda d: f"{abs(d)} a la {'derecha' if d > 0 else 'izquierda'}"
    print(f"f({b_*x + c_}) = f({b_}·(x - ({h})))   ->   h = {h} ≈ {cifras(h, 5)}")
    print(f"Orden A: escalar x por 1/{b_} y después mover {lado(h)}.")
    print(f"Orden B: mover {lado(-c_)} y después escalar x por 1/{b_}.")
    bien = f.subs(x, b_ * x + c_)
    mal = f.subs(x, b_ * (x + c_))
    print("¿mover c en lugar de h da lo mismo?", "sí" if sp.simplify(bien - mal) == 0 else "no")


widgets.interact(calculadora_orden,
    base=widgets.Dropdown(options=list(BASES), value="x²", description="base"),
    b=widgets.FloatText(value=2.0, description="b"), c=widgets.FloatText(value=-6.0, description="c"));

# %%
# Casos de prueba de la sección 4 (resultado conocido)
PRUEBAS_4 = [
    ("(2x - 6): h = 3", lambda: factorizar_interior(2, -6) == 3),
    ("(4x + 8): h = -2", lambda: factorizar_interior(4, 8) == -2),
    ("sen(2x + π/3): h = -π/6", lambda: factorizar_interior(2, sp.pi/3) == -sp.pi/6),
    ("10 sen(120πt - π/3): h = 1/360 s", lambda: factorizar_interior(120*sp.pi, -sp.pi/3) == sp.Rational(1, 360)),
    ("(2x - 6)² y (2(x - 6))² son distintas", lambda: sp.simplify((2*x - 6)**2 - (2*(x - 6))**2) != 0),
    ("2 sen x + 3 y 2(sen x + 3) tienen rangos [1,5] y [4,8]",
     lambda: function_range(2*sp.sin(x) + 3, x, sp.S.Reals) == sp.Interval(1, 5)
             and function_range(2*(sp.sin(x) + 3), x, sp.S.Reals) == sp.Interval(4, 8)),
    ("orden B: mover π/3 y luego comprimir da sen(2x + π/3)",
     lambda: sp.simplify(sp.sin(x + sp.pi/3).subs(x, 2*x) - sp.sin(2*x + sp.pi/3)) == 0),
    ("b = 0 se rechaza", lambda: _rechaza(lambda: factorizar_interior(0, 1))),
]
for nombre, prueba in PRUEBAS_4:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 4).** Escribe a mano $(3x+12)^2$ en la forma $f\big(b(x-h)\big)$ con $f(x)=x^2$. Da $b$ y $h$ y compáralos con la calculadora ($b=3$, $c=12$).

# %%
mi_b, mi_h = None, None      # escribe dos números, por ejemplo: 1, 0

ref = factorizar_interior(3, 12)
print("La calculadora da: b = 3, h =", ref)
print("falta tu cálculo" if mi_h is None else ("coincide" if cerca(mi_h, float(ref)) and mi_b == 3 else "NO coincide: h = -c/b"))
