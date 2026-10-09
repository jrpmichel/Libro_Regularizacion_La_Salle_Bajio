# ID: INT-U5-NB01
# Notebook: int/u5_tecnicas.ipynb · sección 5.1 sustitución
# Repositorio: int/u5_tecnicas/01_sustitucion.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Sustitución
#
# $$\int f\big(g(x)\big)\,g'(x)\,dx=\int f(u)\,du=F\big(g(x)\big)+C,\qquad u=g(x),\quad du=g'(x)\,dx.$$
#
# En una integral definida, los límites cambian con la variable: $\int_a^b f\big(g(x)\big)g'(x)\,dx=\int_{g(a)}^{g(b)}f(u)\,du$.
#
# Primero, qué técnica elige `sympy` para cada integral de la semilla. La función `regla_sympy` lee la estructura de pasos de `integral_steps`.

# %%
semilla = [2 * x * sp.cos(x**2), x * sp.cos(x), 1 / (x**2 - 1), sp.sin(x)**2]
for f in semilla:
    print(f"∫ {f} dx  →  {regla_sympy(f)}")
print("∫_1^∞ dx/x² =", sp.integrate(1 / x**2, (x, 1, sp.oo)), " (impropia, sección 5)")

# %% [markdown]
# Un ejemplo con cambio de límites: $\int_0^1\frac{x}{1+x^2}\,dx$ con $u=1+x^2$, que va de $1$ a $2$.

# %%
directa = sp.integrate(x / (1 + x**2), (x, 0, 1))
en_u = sp.Rational(1, 2) * sp.integrate(1 / U, (U, 1, 2))
print("en x:", directa, "  en u:", en_u, "  ≈", cifras(directa, 4))

# %% [markdown]
# **Calculadora de sustitución.** Escribe $f(x)$ y tu propuesta de $u(x)$. La calculadora calcula $du$ e intenta escribir $f(x)\,dx$ como $h(u)\,du$. Si queda una $x$, intenta despejarla de $u=g(x)$ y te avisa si no lo logra. Si la sustitución sirve, integra en $u$, regresa a $x$ con la constante $C$, deriva el resultado y lo compara con la antiderivada de `sympy`.

# %%
def _nunca_negativa(g):
    """True si g(x) >= 0 en todos los puntos de prueba donde es real (como x² + 1 o √x)."""
    if g.is_nonnegative:
        return True
    valores = []
    for c in _PUNTOS:
        try:
            valores.append(valor_real(g, c))
        except (ValueError, TypeError, ZeroDivisionError, OverflowError):
            continue
    return len(valores) >= 3 and min(valores) >= 0


def _integra_en_u(hu):
    """∫ h(u) du con las reglas "a mano" de sympy (manualintegrate): rápidas y sin funciones especiales."""
    H = con_limite(lambda: _mi.manualintegrate(sp.expand(hu), U), 8)   # expandido: u^(5/2) en vez de (u - 1)² √u
    if H.has(sp.Integral):
        raise ValueError(f"la integral en u, ∫ {hu} du, no está en la tabla: esa sustitución no la simplifica.")
    H = _a_logaritmos(H)
    raras = [a for a in H.atoms(sp.Function) if not isinstance(a, _PERMITIDAS)]
    if H.has(sp.Piecewise) or raras:           # caso raro: se usa la antiderivada general, escrita en u
        H = antiderivada(hu.subs(U, x)).subs(x, U)
    return H


def sustituir(f, g):
    """∫ f(x) dx con u = g(x). Devuelve du/dx, h(u) con f dx = h du, H(u), F(x) y la lista de despejes usados.
    Si al despejar x hay varias ramas (x = ±√u), se queda con la que cumple F' = f.
    ValueError si queda una x suelta, si la integral en u no está en la tabla o si ninguna rama sirve."""
    if not g.has(x):
        raise ValueError("u debe depender de x.")
    if sp.simplify(g - x) == 0:
        raise ValueError("u = x no cambia la integral: elige como u una parte del integrando "
                         "(lo de adentro de una raíz, de un coseno o de una exponencial).")
    dg = sp.diff(g, x)
    if sp.simplify(dg) == 0:
        raise ValueError("u es constante: su derivada es cero.")
    h = sp.simplify(f / dg)                    # f dx = (f / g') du
    candidatos = []
    hu = sp.simplify(h.subs(g, U))
    if not hu.has(x):
        candidatos.append((hu, None))
    else:                                      # queda una x: se intenta despejar de u = g(x)
        try:
            inversas = con_limite(lambda: sp.solve(sp.Eq(U, g), x), 8)
        except (ValueError, NotImplementedError):
            inversas = []
        for inv in inversas:
            if inv.has(sp.I):                  # raíces complejas de u = x³: no sirven en los reales
                continue
            cand = sp.simplify(h.subs(x, inv))
            if cand.has(x):
                continue
            previo = [k for k, (hc, _) in enumerate(candidatos) if sp.simplify(hc - cand) == 0]
            if previo:                         # otra rama da el mismo h(u): x solo aparece al cuadrado
                candidatos[previo[0]][1].append(inv)
            else:
                candidatos.append((cand, [inv]))
    if not candidatos:
        raise ValueError("queda una x que no se puede escribir con u: esa sustitución no sirve.")
    positiva, error = _nunca_negativa(g), None
    for hu, despeje in candidatos:
        if positiva:                           # u >= 0: |u| = u
            hu = hu.subs(sp.Abs(U), U)
        try:
            H = _integra_en_u(hu)
        except ValueError as err:
            error = err; continue
        F = con_valor_absoluto(sp.simplify(H.subs(U, g)))
        if misma_funcion(derivada_legible(F), sin_abs_en_log(f)):
            return {"du": dg, "hu": hu, "H": H, "F": F, "despeje": despeje}
    if error is not None and len(candidatos) == 1:
        raise error
    raise ValueError("al despejar x de u = g(x) hay que elegir un signo (x = ±...) y ninguna elección sirve para "
                     "todo el dominio: esa sustitución no funciona aquí.")


def calculadora_sustitucion(f_txt, u_txt):
    try:
        f, g = parsear(f_txt), parsear(u_txt)
        r = sustituir(f, g)
    except ValueError as err:
        print("No sirve o hay un error en la entrada:", err); return
    print(f"u = {g}     du = ({r['du']}) dx")
    if r["despeje"]:
        ramas = " o x = ".join(str(d) for d in r["despeje"])
        print(f"quedaba una x; de u = {g} se despejó x = {ramas}"
              + ("  (las dos ramas dan el mismo integrando en u)" if len(r["despeje"]) > 1 else ""))
    print(f"∫ f(x) dx = ∫ {r['hu']} du = {r['H']} + C")
    print(f"regreso a x: {mostrar_antiderivada(r['F'])}")
    ok, cte = compara_con_sympy(r["F"], f)
    print("derivada de tu resultado = f(x):", "sí" if ok else "NO")
    if texto_constante(cte):
        print(texto_constante(cte))
    print("regla que eligió sympy:", regla_sympy(f))


widgets.interact(calculadora_sustitucion,
    f_txt=widgets.Text(value="x*cos(x^2)", description="f(x) =", continuous_update=False),
    u_txt=widgets.Text(value="x^2", description="u(x) =", continuous_update=False));

# %%
# Casos de prueba de la sección 1 (resultado conocido)
PRUEBAS_1 = [
    ("x cos x² con u = x²: ½ sen x²", lambda: cerca(valor_real(sustituir(x * sp.cos(x**2), x**2)["F"] - sp.sin(x**2) / 2, 0.7), 0)),
    ("x√(x+1) con u = x+1: despeja x = u - 1",
     lambda: (lambda r: bool(r["despeje"]) and compara_con_sympy(r["F"], x * sp.sqrt(x + 1))[0])(
         sustituir(x * sp.sqrt(x + 1), x + 1))),
    ("cos x³ con u = x³: no sirve (integral en u fuera de la tabla)", lambda: _rechaza(lambda: sustituir(sp.cos(x**3), x**3))),
    ("x cos x³ con u = x³: no sirve", lambda: _rechaza(lambda: sustituir(x * sp.cos(x**3), x**3))),
    ("tan x con u = cos x: -ln|cos x|",
     lambda: compara_con_sympy(sustituir(sp.tan(x), sp.cos(x))["F"], sp.tan(x))[0]
     and sustituir(sp.tan(x), sp.cos(x))["F"].has(sp.Abs)),
    ("2x² con u = x²: no sirve (x = ±√u, ninguna rama vale en todo el dominio)",
     lambda: _rechaza(lambda: sustituir(2 * x**2, x**2))),
    ("u = x se rechaza", lambda: _rechaza(lambda: sustituir(x * sp.cos(x), x))),
]
for nombre, prueba in PRUEBAS_1:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 1).** Calcula a mano $\int x^2e^{x^3}\,dx$ y escribe tu antiderivada como texto, sin la $C$ (por ejemplo `"x**2/2"`).

# %%
mi_respuesta = None     # escribe tu antiderivada entre comillas

if mi_respuesta is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    f = x**2 * sp.exp(x**3)
    try:
        ok, cte = compara_con_sympy(parsear(mi_respuesta), f)
    except ValueError as err:
        print("Revisa tu respuesta:", err)
    else:
        print("La calculadora da:", mostrar_antiderivada(sustituir(f, x**3)["F"]))
        print("tu respuesta es correcta" + (f" (difiere en la constante {cte})" if cte not in (None, 0) else "")
              if ok else "NO coincide: con u = x³, du = 3x² dx, así que x² dx = du/3")
