# ID: DIF-U2-NB04
# Notebook: dif/u2_limite.ipynb · sección 2.4 álgebra de límites y la forma 0/0
# Repositorio: dif/u2_limite/04_algebra.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 4. Álgebra de límites y la forma indeterminada $0/0$
#
# La estrategia del libro, en tres pasos:
#
# 1. **Sustituye** $x=a$. Si sale un número y la función es de las que se pueden evaluar directamente (polinomios, cocientes con denominador distinto de cero, raíces, exponenciales, logaritmos, seno y coseno dentro de su dominio), ese número es el límite.
# 2. Si sale $k/0$ con $k\neq0$, no hay límite finito: revisa los límites laterales, que suelen ser $\pm\infty$.
# 3. Si sale $0/0$, **todavía no sabes nada**. Factoriza y cancela, o multiplica por el conjugado si hay una raíz, y vuelve a sustituir en la expresión simplificada.
#
# La calculadora sigue esos pasos y muestra cada uno. La simplificación que propone es una ayuda: compárala con la tuya.

# %%
def _racionaliza(f):
    """Si el numerador es A + B con una raíz, multiplica arriba y abajo por A - B."""
    num, den = sp.fraction(sp.together(f))
    terminos = sp.Add.make_args(num)
    tiene_raiz = lambda t: any(p.exp == sp.Rational(1, 2) for p in t.atoms(sp.Pow))
    if len(terminos) == 2 and any(tiene_raiz(t) for t in terminos):
        R, C = terminos if tiene_raiz(terminos[0]) else terminos[::-1]   # R: término con raíz
        conj = R - C                                 # (R + C)(R - C) = R² - C²: la raíz desaparece
        return sp.expand(num * conj), den, conj
    return None


def estrategia(texto, a):
    """Devuelve (lista de pasos en texto, límite izquierdo, límite derecho)."""
    f = parsear(texto)
    a = numero(a)
    pasos = []
    valor = f.subs(x, a)
    num, den = sp.fraction(sp.together(f))
    n_a, d_a = sp.simplify(num.subs(x, a)), sp.simplify(den.subs(x, a))
    if valor.is_finite and d_a != 0:
        pasos.append(f"1. Sustitución directa: f({a}) = {a_texto(valor)}. Es un número: ese es el límite.")
    elif n_a != 0 and d_a == 0:
        pasos.append(f"1. Sustitución directa: {n_a}/0. No hay límite finito; mira los laterales (abajo).")
    elif n_a == 0 and d_a == 0:
        pasos.append("1. Sustitución directa: 0/0. Forma indeterminada: todavía no hay respuesta.")
        g = sp.cancel(f)
        if g != f and g.subs(x, a).is_finite:
            pasos.append(f"2. Factorizo: ({sp.factor(num)}) / ({sp.factor(den)}). Cancelo el factor que se anula "
                         f"(se vale porque x ≠ {a}): {g}.")
        else:
            r = _racionaliza(f)
            if r is not None:
                nuevo_num, den0, conj = r
                g = sp.cancel(nuevo_num / den0) / conj
                pasos.append(f"2. Multiplico arriba y abajo por el conjugado {conj}: el numerador queda {nuevo_num}; "
                             f"cancelo y obtengo {g}.")
            else:
                pasos.append("2. No encontré una simplificación algebraica directa; usa una tabla (sección 1) "
                             "o los límites trigonométricos (sección 5).")
        if g.subs(x, a).is_finite:
            pasos.append(f"3. Sustituyo en la expresión simplificada: {a_texto(g.subs(x, a))}.")
    else:
        pasos.append(f"1. f({a}) no está definida ({valor}). Revisa el dominio.")
    li, ld = sp.limit(f, x, a, "-"), sp.limit(f, x, a, "+")
    return pasos, li, ld


def calculadora_algebra(texto, a):
    try:
        pasos, li, ld = estrategia(texto, a)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"lím f(x) con x -> {a},   f(x) = {texto}")
    for p in pasos:
        print("  ", p)
    print(f"   Comprobación con sympy: izquierda {a_texto(li)}, derecha {a_texto(ld)}.")


widgets.interact(calculadora_algebra,
    texto=widgets.Text(value="(x**3 - 1)/(x - 1)", description="f(x) ="),
    a=widgets.Text(value="1", description="a ="));

# %% [markdown]
# **Una aplicación con $0/0$: la diferencia media logarítmica de temperatura.** En un intercambiador de calor se usa $\Delta T_{\mathrm{ml}}=\dfrac{\Delta T_1-\Delta T_2}{\ln(\Delta T_1/\Delta T_2)}$. Cuando las dos diferencias son iguales, la fórmula da $0/0$. La tabla muestra hacia dónde va con $\Delta T_1=20$ °C. El cálculo exacto llega en la Unidad 6.

# %%
dTml = lambda d1, d2: (d1 - d2) / math.log(d1 / d2)
for d2 in (19.0, 19.9, 19.99, 20.01, 20.1, 21.0):
    print(f"ΔT2 = {d2:<6} °C   ΔT_ml = {dTml(20.0, d2):.5f} °C")

# %%
# Casos de prueba de la sección 4 (resultado conocido)
def _lim(texto, a):
    _, li, ld = estrategia(texto, a)
    return li if li == ld else None

PRUEBAS_4 = [
    ("(x³-1)/(x-1) en 1 da 3", lambda: _lim("(x**3 - 1)/(x - 1)", "1") == 3),
    ("(x²-9)/(x-3) en 3 da 6 (0/0 no significa 'no existe')", lambda: _lim("(x**2 - 9)/(x - 3)", "3") == 6),
    ("(sqrt(x+1)-1)/x en 0 da 1/2, con conjugado",
     lambda: _lim("(sqrt(x + 1) - 1)/x", "0") == sp.Rational(1, 2)
             and any("conjugado" in p for p in estrategia("(sqrt(x + 1) - 1)/x", "0")[0])),
    ("(sqrt(x+3)-2)/(x²-1) en 1 da 1/8 (problema DIF-U2-03)", lambda: _lim("(sqrt(x + 3) - 2)/(x**2 - 1)", "1") == sp.Rational(1, 8)),
    ("(x+1)/(x-3) en 3: 4/0, laterales -∞ y +∞",
     lambda: estrategia("(x + 1)/(x - 3)", "3")[1:] == (-sp.oo, sp.oo)),
    ("x²+1 en 2: sustitución directa da 5", lambda: "Sustitución directa: f(2) = 5" in estrategia("x**2 + 1", "2")[0][0]),
    ("ΔT_ml con ΔT2 = 19.99 está a menos de 0.01 °C de 20", lambda: abs(dTml(20.0, 19.99) - 20) < 0.01),
]
for nombre, prueba in PRUEBAS_4:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 4).** Calcula a mano $\displaystyle\lim_{x\to3}\frac{x^2-x-6}{x-3}$: sustituye, factoriza el numerador, cancela y vuelve a sustituir. Escribe tu resultado y compáralo con la calculadora. Anota qué factor cancelaste y por qué se vale cancelarlo.

# %%
mi_limite = None        # por ejemplo: 5

ref = _lim("(x**2 - x - 6)/(x - 3)", "3")
print("La calculadora da:", ref)
print("falta tu cálculo" if mi_limite is None else ("coincide" if cerca(mi_limite, float(ref)) else "NO coincide: revisa la factorización"))
