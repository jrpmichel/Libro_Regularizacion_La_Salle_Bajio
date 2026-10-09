# ID: INT-U5-NB02
# Notebook: int/u5_tecnicas.ipynb · sección 5.2 integración por partes
# Repositorio: int/u5_tecnicas/02_partes.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Integración por partes
#
# $$\int u\,dv=uv-\int v\,du,\qquad\int_a^b u\,v'\,dx=\Big[uv\Big]_a^b-\int_a^b v\,u'\,dx.$$
#
# Sale de la regla del producto, $(uv)'=u'v+uv'$. La figura del libro lo muestra con áreas: para $y=\ln x$ en $[1,e]$, el área bajo la curva más el área a su izquierda llenan el rectángulo $e\times1$.

# %%
bajo = sp.integrate(sp.log(x), (x, 1, sp.E))
izquierda = sp.integrate(sp.exp(U), (U, 0, 1))
print("bajo ln x:", bajo, "  a la izquierda:", izquierda, "  suma:", sp.simplify(bajo + izquierda), "= e")

# %% [markdown]
# **Calculadora de partes.** Escribe $f(x)$ y tu elección de $u$; el resto es $dv=\frac{f}{u}\,dx$. La calculadora muestra $du$, $v$ y la fórmula $\int f\,dx=uv-\int v\,du$ con la integral nueva, y te avisa si esa integral es más complicada que la original: si la contiene, si sube el grado del polinomio, si aparece una función trascendente que no estaba o, con lo demás igual, si tiene más operaciones. Si la integral nueva es un polinomio por una exponencial o una función trigonométrica, aplica partes otra vez, hasta tres veces.

# %%
_CLASE = {sp.exp: "exp", sp.log: "ln", sp.sin: "trig", sp.cos: "trig", sp.tan: "trig", sp.sec: "trig",
          sp.atan: "arco", sp.asin: "arco", sp.acos: "arco"}


def _tipos(g):
    """Clases de funciones trascendentes en g: exp, ln, trig (sen, cos y tan cuentan igual) y arco."""
    return {_CLASE.get(type(a), type(a).__name__) for a in g.atoms(sp.Function)}


def _grado_polinomial(g):
    """Grado del factor polinomial de g (0 si no tiene)."""
    poli = sp.Mul(*[a for a in sp.Mul.make_args(sp.factor(g)) if a.has(x) and a.is_polynomial(x)])
    return sp.degree(poli, x) if poli.has(x) else 0


def _empeora(nueva, f):
    """La integral nueva es peor si contiene a la original, si aparece una función trascendente que no estaba,
    si sube el grado del polinomio
    o, con las mismas funciones y el mismo grado, si tiene más operaciones. Si desaparece una función
    trascendente (ln x o arctan x se derivan y se van), la integral nueva es mejor."""
    if any(not sp.simplify(t / f).has(x) for t in sp.Add.make_args(sp.expand(nueva))):
        return True                            # la integral nueva contiene a la original: no se avanzó
    tn, tf = _tipos(nueva), _tipos(f)
    if tn < tf:
        return False
    if not tn <= tf:
        return True
    gn, gf = _grado_polinomial(nueva), _grado_polinomial(f)
    return gn > gf or (gn == gf and sp.count_ops(nueva) > sp.count_ops(f) + 1)


def partes(f, u_expr):
    """Un paso: u = u_expr, dv = f/u. Devuelve u, du, dv, v, la integral nueva v du y si la elección empeoró la integral."""
    if not u_expr.has(x):
        raise ValueError("u debe depender de x.")
    dv = sp.simplify(f / u_expr)
    try:
        v = antiderivada(dv)
    except ValueError:
        raise ValueError(f"con u = {u_expr} queda dv = ({dv}) dx, y ese dv no se sabe integrar: elige otra u.") from None
    du = sp.diff(u_expr, x)
    nueva = sp.simplify(v * du)
    return {"f": f, "u": u_expr, "du": du, "dv": dv, "v": v, "nueva": nueva, "empeora": _empeora(nueva, f)}


def _parte_polinomial(g):
    """Si g es un polinomio de grado >= 1 por una exponencial o una función trigonométrica, el polinomio
    (la siguiente u); si no, None."""
    fac = sp.Mul.make_args(sp.factor(g))
    poli = sp.Mul(*[a for a in fac if a.has(x) and a.is_polynomial(x)])
    resto = sp.Mul(*[a for a in fac if not (a.has(x) and a.is_polynomial(x))])
    sin_cocientes = not any(isinstance(a, sp.Pow) and a.exp.is_negative for a in sp.Mul.make_args(resto))
    if poli.has(x) and _tipos(resto) and _tipos(resto) <= {"exp", "trig"} and sin_cocientes:
        return poli
    return None


def por_partes(f, u_expr, max_pasos=3):
    """Aplica partes hasta max_pasos veces. Devuelve la lista de pasos y la antiderivada F (con C = 0).
    El primer paso usa la u del alumno; los siguientes, el factor polinomial de la integral nueva."""
    pasos, signo, acumulado, actual, u_act = [], 1, 0, f, u_expr
    for k in range(max_pasos):
        try:
            p = partes(actual, u_act)
        except ValueError:
            if k == 0:
                raise
            break                              # el paso siguiente no se puede: se integra lo que queda
        pasos.append(p)
        acumulado += signo * p["u"] * p["v"]
        signo, actual = -signo, p["nueva"]     # F = uv - ∫ v du: la integral nueva entra con el signo cambiado
        if p["empeora"]:
            break
        u_act = _parte_polinomial(actual)
        if u_act is None:
            break
    try:
        resto = antiderivada(actual)
    except ValueError:
        raise ValueError(f"la integral que queda, ∫ {actual} dx, no tiene una antiderivada elemental; "
                         "prueba otra elección de u.") from None
    return pasos, sin_desfase(sp.simplify(acumulado + signo * resto))


def calculadora_partes(f_txt, u_txt):
    try:
        f, u_expr = parsear(f_txt), parsear(u_txt)
        pasos, F = por_partes(f, u_expr)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    for k, p in enumerate(pasos, 1):
        print(f"paso {k}: u = {p['u']}, dv = ({p['dv']}) dx  →  du = ({p['du']}) dx, v = {p['v']}")
        print(f"        ∫ {p['f']} dx = {p['u'] * p['v']} - ∫ {p['nueva']} dx")
        if p["empeora"]:
            print("        AVISO: la integral nueva es más complicada que la original; prueba otra elección de u")
    print("resultado:", mostrar_antiderivada(F))
    ok, cte = compara_con_sympy(F, f)
    print("derivada del resultado = f(x):", "sí" if ok else "NO")
    if texto_constante(cte):
        print(texto_constante(cte))
    print("regla que eligió sympy:", regla_sympy(f))


widgets.interact(calculadora_partes,
    f_txt=widgets.Text(value="x*exp(x)", description="f(x) =", continuous_update=False),
    u_txt=widgets.Text(value="x", description="u =", continuous_update=False));

# %%
# Casos de prueba de la sección 2 (resultado conocido)
PRUEBAS_2 = [
    ("x e^x con u = x: (x - 1)e^x", lambda: sp.simplify(por_partes(x * sp.exp(x), x)[1] - (x - 1) * sp.exp(x)) == 0),
    ("x e^x con u = e^x: la integral nueva empeora", lambda: partes(x * sp.exp(x), sp.exp(x))["empeora"]),
    ("ln x con u = ln x: x ln x - x", lambda: misma_funcion(por_partes(sp.log(x), sp.log(x))[1], x * sp.log(x) - x)),
    ("x² sen x con u = x²: dos pasos, -x² cos x + 2x sen x + 2 cos x",
     lambda: (lambda r: len(r[0]) == 2 and sp.simplify(r[1] - (-x**2 * sp.cos(x) + 2 * x * sp.sin(x) + 2 * sp.cos(x))) == 0)(
         por_partes(x**2 * sp.sin(x), x**2))),
    ("x ln x con u = ln x: la integral nueva (x/2) no empeora", lambda: not partes(x * sp.log(x), sp.log(x))["empeora"]),
    ("x e^x con u = x²: dv = e^x/x no se integra", lambda: _rechaza(lambda: partes(x * sp.exp(x), x**2))),
]
for nombre, prueba in PRUEBAS_2:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 2).** Calcula a mano $\int x\cos x\,dx$ y escribe tu antiderivada como texto, sin la $C$.

# %%
mi_respuesta = None     # escribe tu antiderivada entre comillas

if mi_respuesta is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    f = x * sp.cos(x)
    try:
        ok, cte = compara_con_sympy(parsear(mi_respuesta), f)
    except ValueError as err:
        print("Revisa tu respuesta:", err)
    else:
        print("La calculadora da:", mostrar_antiderivada(por_partes(f, x)[1]))
        print("tu respuesta es correcta" + (f" (difiere en la constante {cte})" if cte not in (None, 0) else "")
              if ok else "NO coincide: con u = x y dv = cos x dx, v = sen x; revisa el signo de ∫ sen x dx")
