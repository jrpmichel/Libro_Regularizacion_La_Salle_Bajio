# ID: INT-U5-NB03
# Notebook: int/u5_tecnicas.ipynb · sección 5.3 fracciones parciales
# Repositorio: int/u5_tecnicas/03_fracciones.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Fracciones parciales
#
# Una función racional propia $\frac{P}{Q}$ se escribe como suma de fracciones simples: $\frac{A}{x-r}$ por cada factor lineal, $\frac{A_1}{x-r}+\dots+\frac{A_k}{(x-r)^k}$ por cada factor lineal repetido y $\frac{Bx+C}{x^2+bx+c}$ por cada factor cuadrático irreducible. Si no es propia, primero se divide.
#
# `sympy.apart` hace la descomposición. Aquí se comprueba con el ejemplo del libro, $\frac{2}{x^2-1}$.

# %%
f = 2 / (x**2 - 1)
piezas = sp.apart(f, x)
print("piezas:", piezas, "   suman f:", sp.simplify(piezas - f) == 0)
print("∫_2^3 =", sp.simplify(sp.integrate(f, (x, 2, 3))), "≈", cifras(sp.integrate(f, (x, 2, 3)), 6))

# %% [markdown]
# **Calculadora de fracciones parciales.** Escribe una función racional. La calculadora revisa los grados (y divide si hace falta), factoriza el denominador, muestra la descomposición, comprueba que las piezas suman la original e integra término por término. Si hay un factor cuadrático irreducible, te muestra cómo completar el cuadrado. Los denominadores con factores de grado mayor que 2 que no se factorizan con coeficientes racionales quedan fuera del curso y la calculadora los rechaza.

# %%
def _integra_pieza(t):
    """∫ c·g dx = c ∫ g dx: la constante se separa, para que salga ¼ ln|x - 2| y no ¼ ln|4x - 8|."""
    c, g = t.as_independent(x, as_Add=False)
    return c * antiderivada(g)


def fracciones(f):
    """Descomposición e integral de una función racional. ValueError si f no es un cociente de polinomios
    o si el denominador tiene un factor de grado mayor que 2 sin factorizar."""
    num, den = sp.fraction(sp.cancel(sp.together(f)))
    if not (num.is_polynomial(x) and den.is_polynomial(x)):
        raise ValueError("f debe ser un cociente de polinomios en x.")
    if not den.has(x):
        raise ValueError(f"al simplificar, f = {sp.expand(num / den)} es un polinomio (no queda x en el denominador): "
                         "intégralo con la tabla de la Unidad 3.")
    constante, factores = sp.factor_list(den)
    grandes = [g for g, _ in factores if sp.degree(g, x) > 2]
    if grandes:
        raise ValueError(f"el factor {grandes[0]} del denominador tiene grado {sp.degree(grandes[0], x)} y no se "
                         "factoriza con coeficientes racionales: ese caso está fuera del curso.")
    cociente, resto = sp.div(sp.Poly(num, x), sp.Poly(den, x))
    propia = sp.degree(num, x) < sp.degree(den, x)
    cuadraticos = [g for g, _ in factores if sp.degree(g, x) == 2 and sp.discriminant(g, x) < 0]
    irracionales = [g for g, _ in factores if sp.degree(g, x) == 2 and sp.discriminant(g, x) > 0]
    piezas = sp.apart(resto.as_expr() / den, x)
    terminos = [t for t in [cociente.as_expr(), *sp.Add.make_args(piezas)] if t != 0]
    F = sp.Add(*[_integra_pieza(t) for t in terminos])
    return {"propia": propia, "cociente": cociente.as_expr(), "constante": constante, "factores": factores,
            "cuadraticos": cuadraticos,
            "irracionales": irracionales, "piezas": piezas, "F": F,
            "suman": sp.simplify(cociente.as_expr() + piezas - f) == 0}


def completar_cuadrado(g):
    """a x² + b x + c -> a·[(x + h)² + k] con h = b/(2a) y k = c/a - h², escrito para leerse."""
    a, b, c = sp.Poly(g, x).all_coeffs()
    h, k = sp.nsimplify(b / (2 * a)), sp.nsimplify(c / a - (b / (2 * a))**2)
    cuadrado = "x²" if h == 0 else f"(x {'+' if h > 0 else '-'} {abs(h)})²"
    texto = f"{cuadrado} {'+' if k >= 0 else '-'} {abs(k)}"
    return texto if a == 1 else f"{a}·[{texto}]"


def calculadora_fracciones(f_txt):
    try:
        f = parsear(f_txt)
        r = fracciones(f)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    if not r["propia"]:
        print(f"no es propia: la división da {r['cociente']} más una fracción propia")
    piezas_den = [f"({g})^{k}" if k > 1 else f"({g})" for g, k in r["factores"]]
    print("denominador factorizado:", " · ".join(([str(r["constante"])] if r["constante"] != 1 else []) + piezas_den))
    for g in r["cuadraticos"]:
        print(f"factor cuadrático irreducible {g} = {completar_cuadrado(g)}  →  da un logaritmo y un arco tangente")
    for g in r["irracionales"]:
        raices = ", ".join(str(v) for v in sp.solve(g, x))
        print(f"AVISO: {g} tiene raíces reales irracionales ({raices}); con coeficientes racionales la descomposición "
              "lo deja junto, y su integral da logaritmos")
    print("descomposición:", r["piezas"] if r["cociente"] == 0 else f"{r['cociente']} + ({r['piezas']})")
    print("las piezas suman la función original:", "sí" if r["suman"] else "NO")
    print("∫ f dx =", mostrar_antiderivada(r["F"]))
    ok, cte = compara_con_sympy(r["F"], f)
    print("derivada del resultado = f(x):", "sí" if ok else "NO")
    if texto_constante(cte):
        print(texto_constante(cte))
    print("regla que eligió sympy:", regla_sympy(f))


widgets.interact(calculadora_fracciones,
    f_txt=widgets.Text(value="(x+5)/(x^2+x-2)", description="f(x) =", continuous_update=False));

# %%
# Casos de prueba de la sección 3 (resultado conocido)
PRUEBAS_3 = [
    ("(x+5)/(x²+x-2) = 2/(x-1) - 1/(x+2)",
     lambda: sp.simplify(fracciones((x + 5) / (x**2 + x - 2))["piezas"] - (2 / (x - 1) - 1 / (x + 2))) == 0),
    ("1/(x(x-1)²) = 1/x - 1/(x-1) + 1/(x-1)²",
     lambda: sp.simplify(fracciones(1 / (x * (x - 1)**2))["piezas"] - (1 / x - 1 / (x - 1) + 1 / (x - 1)**2)) == 0),
    ("(2x+3)/(x²+2x+5): factor irreducible y F' = f",
     lambda: (lambda r: len(r["cuadraticos"]) == 1 and compara_con_sympy(r["F"], (2 * x + 3) / (x**2 + 2 * x + 5))[0])(
         fracciones((2 * x + 3) / (x**2 + 2 * x + 5)))),
    ("x³/(x²-1): no es propia, cociente x", lambda: (lambda r: not r["propia"] and r["cociente"] == x)(fracciones(x**3 / (x**2 - 1)))),
    ("sen x se rechaza (no es racional)", lambda: _rechaza(lambda: fracciones(sp.sin(x)))),
    ("1/(x³ - 2) se rechaza (factor de grado 3)", lambda: _rechaza(lambda: fracciones(1 / (x**3 - 2)))),
    ("1/(x² - 4): ¼ ln|x - 2| - ¼ ln|x + 2|, sin constantes sueltas",
     lambda: sp.simplify(fracciones(1 / (x**2 - 4))["F"]
                         - (sp.log(sp.Abs(x - 2)) / 4 - sp.log(sp.Abs(x + 2)) / 4)) == 0),
]
for nombre, prueba in PRUEBAS_3:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 3).** Calcula a mano $\int\frac{dx}{x^2-4}$ y escribe tu antiderivada como texto, sin la $C$.

# %%
mi_respuesta = None     # escribe tu antiderivada entre comillas

if mi_respuesta is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    f = 1 / (x**2 - 4)
    try:
        ok, cte = compara_con_sympy(parsear(mi_respuesta), f)
    except ValueError as err:
        print("Revisa tu respuesta:", err)
    else:
        print("La calculadora da:", mostrar_antiderivada(fracciones(f)["F"]))
        print("tu respuesta es correcta" + (f" (difiere en la constante {cte})" if cte not in (None, 0) else "")
              if ok else "NO coincide: 1/(x² - 4) = (1/4)/(x - 2) - (1/4)/(x + 2)")
