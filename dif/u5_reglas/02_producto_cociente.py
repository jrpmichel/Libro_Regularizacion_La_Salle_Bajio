# ID: DIF-U5-NB02
# Notebook: dif/u5_reglas.ipynb · sección 5.2 producto y cociente
# Repositorio: dif/u5_reglas/02_producto_cociente.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Reglas de producto y cociente
#
# $$\big(u\,v\big)'=u'\,v+u\,v',\qquad \left(\frac{u}{v}\right)'=\frac{u'\,v-u\,v'}{v^2}\quad (v\neq0).$$
#
# La derivada de un producto **no** es el producto de las derivadas: en el rectángulo de lados $u$ y $v$, el área cambia por dos tiras, $v\,\Delta u$ y $u\,\Delta v$, y una esquina $\Delta u\,\Delta v$ que desaparece en el límite. En el cociente, el orden del numerador importa: invertirlo cambia el signo.
#
# La calculadora recibe $u(x)$ y $v(x)$, arma $(uv)'$ y $(u/v)'$ con las reglas, compara con `sp.diff` y avisa dónde $v(x)=0$.

# %%
def producto_cociente(u_txt, v_txt):
    """Diccionario con u', v', (uv)' y (u/v)' por las reglas, y los ceros reales de v."""
    U, V = parsear(u_txt), parsear(v_txt)
    if V == 0:
        raise ValueError("v(x) = 0 para todo x: el cociente no está definido.")
    dU, dV = sp.diff(U, x), sp.diff(V, x)
    prod = sp.simplify(dU * V + U * dV)
    coc = sp.simplify((dU * V - U * dV) / V**2)
    ceros = sorted(sp.solveset(V, x, domain=sp.S.Reals), key=float) if V.is_polynomial(x) else None
    return {"u'": dU, "v'": dV, "(uv)'": prod, "(u/v)'": coc, "u'v'": sp.simplify(dU * dV), "ceros_v": ceros,
            "ok": iguales(prod, sp.diff(U * V, x)) and iguales(coc, sp.diff(U / V, x))}


def calculadora_producto(u_txt, v_txt):
    try:
        r = producto_cociente(u_txt, v_txt)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print("u' =", r["u'"], "     v' =", r["v'"])
    print("(uv)'  = u'v + uv'          =", r["(uv)'"])
    print("(u/v)' = (u'v - uv')/v^2    =", r["(u/v)'"])
    print("u'v' (NO es la derivada de uv) =", r["u'v'"])
    if r["ceros_v"] is None:
        print("v no es polinomio: revisa a mano dónde se anula.")
    elif r["ceros_v"]:
        print("el cociente no está definido donde v = 0: x =", ", ".join(str(c) for c in r["ceros_v"]))
    print("comprobación con sp.diff:", "coinciden" if r["ok"] else "NO coinciden")


widgets.interact(calculadora_producto,
    u_txt=widgets.Text(value="x**2 + 1", description="u(x) ="),
    v_txt=widgets.Text(value="x**3 - x", description="v(x) ="));

# %%
# Casos de prueba de la sección 2 (resultado conocido)
PRUEBAS_2 = [
    ("(x²+1)(x³-x) -> 5x^4 - 1", lambda: iguales(producto_cociente("x**2 + 1", "x**3 - x")["(uv)'"], 5*x**4 - 1)),
    ("en x = 2: u'v + uv' = 79, u'v' = 44",
     lambda: producto_cociente("x**2 + 1", "x**3 - x")["(uv)'"].subs(x, 2) == 79
             and producto_cociente("x**2 + 1", "x**3 - x")["u'v'"].subs(x, 2) == 44),
    ("(x+1)/(x-1) -> -2/(x-1)²", lambda: iguales(producto_cociente("x + 1", "x - 1")["(u/v)'"], -2/(x - 1)**2)),
    ("ceros de v = x - 1: [1]", lambda: producto_cociente("x + 1", "x - 1")["ceros_v"] == [1]),
    ("costo promedio C/q en 100: -0.15",
     lambda: producto_cociente("2000 + 15*x + 0.05*x**2", "x")["(u/v)'"].subs(x, 100) == sp.Rational(-3, 20)),
    ("v = x² + 1 no tiene ceros reales", lambda: producto_cociente("1", "x**2 + 1")["ceros_v"] == []),
    ("v = 0 se rechaza con mensaje", lambda: _rechaza(lambda: producto_cociente("x", "0"))),
]
for nombre, prueba in PRUEBAS_2:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 2).** Deriva a mano $\dfrac{x^2}{x+1}$ con la regla del cociente y escribe el resultado como texto (con el formato de `"(x - 1)/(x + 3)**2"`).

# %%
mi_derivada = None       # escribe un texto entre comillas

contrasta(mi_derivada, sp.diff(x**2 / (x + 1), x),
          pista="¿Escribiste u'v - uv' en ese orden? ¿Elevaste v al cuadrado en el denominador?")
