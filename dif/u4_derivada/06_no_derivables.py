# ID: DIF-U4-NB06
# Notebook: dif/u4_derivada.ipynb · sección 4.6 derivable implica continua; puntos no derivables
# Repositorio: dif/u4_derivada/06_no_derivables.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6. Propiedades: derivable implica continua; puntos no derivables
#
# Las **derivadas laterales** son los límites del cociente incremental con $h\to0^-$ y con $h\to0^+$. $f$ es derivable en $a$ si las dos existen, son finitas e iguales. Si $f$ es derivable en $a$, es continua en $a$; al revés no: $|x|$ es continua en $0$ y no es derivable.
#
# Primero se revisa la continuidad; las filas esquina, cúspide y tangente vertical suponen $f$ continua en $a$.
#
# | Caso | Laterales del cociente | Ejemplo en $x=0$ |
# |---|---|---|
# | discontinuidad | (no hace falta calcularlos) | escalón |
# | esquina | finitos y distintos | $\lvert x\rvert$ |
# | cúspide | $-\infty$ y $+\infty$ (signos opuestos) | $x^{2/3}$ |
# | tangente vertical | los dos $+\infty$ o los dos $-\infty$ | $x^{1/3}$ |
#
# En un extremo del dominio, como $\sqrt{x}$ en $0$, solo existe el cociente de un lado.
#
# La calculadora recibe una función por partes: la regla de la izquierda vale para $x<a$ y la de la derecha para $x\ge a$, así que $f(a)$ sale de la regla de la derecha. Para una sola fórmula escribe la misma en los dos lados. Para $x^{2/3}$ escribe `cbrt(x)**2` y para $x^{1/3}$, `cbrt(x)`.

# %%
def clasificar(izq_txt, der_txt, a):
    """(tipo, lateral izquierda, lateral derecha) de la función por partes en x = a."""
    fi, fd = parsear(izq_txt), parsear(der_txt)
    a = sp.nsimplify(a)
    if not definida_en(fd, a):
        raise ValueError(f"la regla de la derecha no está definida en x = {a}, así que f(a) no existe.")
    fa = fd.subs(x, a)
    if not lado_definido(fi, a, -1):            # extremo izquierdo del dominio: solo cuenta la derecha
        if sp.simplify(sp.limit(fd, x, a, "+") - fa) != 0:
            return "discontinua por la derecha: no es derivable", None, None
        qd = sp.limit((fd.subs(x, a + h) - fa) / h, h, 0, "+")
        if qd in (sp.oo, -sp.oo):
            return "extremo del dominio con tangente vertical (solo existe el lado derecho)", None, qd
        return f"extremo del dominio: solo existe la derivada por la derecha, {qd}", None, qd
    li, ld = sp.limit(fi, x, a, "-"), sp.limit(fd, x, a, "+")
    if not (li.is_finite and ld.is_finite and sp.simplify(li - fa) == 0 and sp.simplify(ld - fa) == 0):
        return "discontinua: no es derivable (derivable implica continua)", None, None
    qi = sp.limit((fi.subs(x, a + h) - fa) / h, h, 0, "-")
    qd = sp.limit((fd.subs(x, a + h) - fa) / h, h, 0, "+")
    if qi.is_finite and qd.is_finite:
        if sp.simplify(qi - qd) == 0:
            return f"derivable, f'(a) = {sp.simplify(qd)}", qi, qd
        return "esquina: laterales finitos y distintos", qi, qd
    if qi in (sp.oo, -sp.oo) and qd in (sp.oo, -sp.oo):
        return ("tangente vertical: los dos laterales con el mismo signo infinito" if qi == qd
                else "cúspide: laterales infinitos de signos opuestos"), qi, qd
    return "no derivable: un lateral finito y otro infinito", qi, qd


def calculadora_laterales(izq_txt, der_txt, a, n_cifras):
    try:
        tipo, qi, qd = clasificar(izq_txt, der_txt, a)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    if qd is not None:
        print("derivada por la izquierda:", lado_texto(qi), "| por la derecha:", a_texto(qd))
    print("veredicto:", tipo)
    fi, fd = sp.lambdify(x, parsear(izq_txt), "numpy"), sp.lambdify(x, parsear(der_txt), "numpy")
    for hh in (0.1, 0.01, 0.001):
        ci = (evaluar(fi, a - hh) - evaluar(fd, a)) / (-hh)
        cd = (evaluar(fd, a + hh) - evaluar(fd, a)) / hh
        ci, cd = (0.0 if abs(c) < 1e-12 else c for c in (ci, cd))   # ruido de punto flotante
        print(f"   h = {hh:<6} cociente izq. {cifras(ci, n_cifras):>10}   der. {cifras(cd, n_cifras):>10}")


widgets.interact(calculadora_laterales,
    izq_txt=widgets.Text(value="abs(x)", description="x < a:"),
    der_txt=widgets.Text(value="abs(x)", description="x ≥ a:"),
    a=widgets.FloatText(value=0.0, description="a"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 6 (resultado conocido)
PRUEBAS_6 = [
    ("|x| en 0: esquina con laterales -1 y 1", lambda: clasificar("abs(x)", "abs(x)", 0)[1:] == (-1, 1)),
    ("x^(2/3) en 0: cúspide", lambda: clasificar("cbrt(x)**2", "cbrt(x)**2", 0)[0].startswith("cúspide")),
    ("x^(1/3) en 0: tangente vertical", lambda: clasificar("cbrt(x)", "cbrt(x)", 0)[0].startswith("tangente vertical")),
    ("x² en 0: derivable con f'(0) = 0", lambda: clasificar("x**2", "x**2", 0)[0] == "derivable, f'(a) = 0"),
    ("x² (x < 1) y 2x - 1 (x ≥ 1): derivable con f'(1) = 2",
     lambda: clasificar("x**2", "2*x - 1", 1)[0] == "derivable, f'(a) = 2"),
    ("escalón 0 / 2 en 0: discontinua", lambda: clasificar("0*x", "2 + 0*x", 0)[0].startswith("discontinua")),
    ("rampa 0 / 400x en 0: esquina (0 y 400)", lambda: clasificar("0*x", "400*x", 0)[1:] == (0, 400)),
    ("√x en 0: extremo del dominio con tangente vertical",
     lambda: clasificar("sqrt(x)", "sqrt(x)", 0)[0].startswith("extremo del dominio con tangente vertical")),
]
for nombre, prueba in PRUEBAS_6:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 6).** Para $f(x)=|x-2|$, calcula a mano las pendientes de las secantes desde $(2,0)$ con $h=0.1$ y con $h=-0.1$. ¿Es derivable en $x=2$? Compara con la calculadora (`abs(x - 2)` en los dos lados y $a=2$).

# %%
mi_izq, mi_der = None, None     # escribe dos números, por ejemplo: 0, 0

ref = clasificar("abs(x - 2)", "abs(x - 2)", 2)
print("La calculadora da:", ref[1], "y", ref[2], "->", ref[0])
print("falta tu cálculo" if mi_izq is None else
      ("coincide" if cerca(mi_izq, float(ref[1])) and cerca(mi_der, float(ref[2])) else "NO coincide: |±0.1|/(±0.1)"))
