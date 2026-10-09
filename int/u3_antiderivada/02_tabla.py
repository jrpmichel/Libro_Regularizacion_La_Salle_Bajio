# ID: INT-U3-NB02
# Notebook: int/u3_antiderivada.ipynb · sección 3.2 tabla de integrales inmediatas y verificador
# Repositorio: int/u3_antiderivada/02_tabla.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 2. Tabla de integrales inmediatas y verificador de antiderivadas
#
# | $f(x)$ | $\int f(x)\,dx$ |
# |---|---|
# | $x^n$, $n\neq-1$ | $\frac{x^{n+1}}{n+1}+C$ |
# | $\frac1x$ | $\ln\lvert x\rvert+C$ |
# | $e^x$ | $e^x+C$ |
# | $a^x$ | $\frac{a^x}{\ln a}+C$ |
# | $\sen x$ | $-\cos x+C$ |
# | $\cos x$ | $\sen x+C$ |
# | $\frac{1}{\cos^2x}$ | $\tan x+C$ |
#
# Regla del argumento lineal: $\int f(ax+b)\,dx=\frac1aF(ax+b)+C$.
#
# El **verificador** recibe $f$ y la antiderivada que tú propones. Deriva tu propuesta y la compara con $f$, primero con `sympy` y, si hace falta, en ocho puntos de $[-3,3]$. Te dice si tu respuesta es correcta, si solo vale en parte del dominio de $f$ (por ejemplo, un `log(x)` sin valor absoluto) o si es incorrecta, y en ese caso cuánto vale $F'-f$. Si tu respuesta es correcta pero se ve distinta de la de `sympy`, te dice en qué constante difieren.

# %%
def verificar(f, F):
    """('correcta' | 'parcial' | 'incorrecta', detalle) para la antiderivada propuesta F de f."""
    dF = derivada_legible(F)
    simbolico = sp.simplify(dF - sin_abs_en_log(f)) == 0
    falla_F, compara, distintos = [], 0, []
    for c in _PUNTOS:
        if not _definida(f, c):
            continue                                  # f no existe ahí: no cuenta
        if not _definida(F, c):
            falla_F.append(c); continue               # f existe pero tu F no
        fc, dc = valor_real(f, c), valor_real(dF, c)
        compara += 1
        if abs(dc - fc) > 1e-9 * max(1, abs(fc)):     # tolerancia relativa
            distintos.append(c)
    if simbolico or (compara >= 3 and not distintos):
        if falla_F:
            lista = ", ".join(cifras(c, 2) for c in falla_F)
            return "parcial", (f"F' = f donde tu F existe, pero tu F no existe en x = {lista}, donde f sí: "
                               "falta un valor absoluto o tu F no cubre todo el dominio de f.")
        G = antiderivada(f)
        dif = sp.simplify(F - G)
        extra = f" Difiere de la de sympy, {G}, en la constante {dif}." if dif != 0 and dif.is_number else ""
        return "correcta", "F' coincide con f." + extra
    detalle = f"F' - f = {sp.simplify(dF - f)}"
    if distintos:
        c = distintos[0]
        detalle += f"; en x = {c} vale {cifras(valor_real(dF, c) - valor_real(f, c), 6)}"
    return "incorrecta", detalle


def _definida(expr, c):
    try:
        valor_real(expr, c); return True
    except (ValueError, TypeError, ZeroDivisionError):
        return False


_TEXTO = {"correcta": "CORRECTA", "parcial": "CORRECTA SOLO EN PARTE DEL DOMINIO", "incorrecta": "INCORRECTA"}


def calculadora_verificador(f_txt, F_txt):
    try:
        f, F = parsear(f_txt), parsear(F_txt)
        veredicto, detalle = verificar(f, F)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"f(x) = {f}\ntu F(x) = {F}\nF'(x) = {derivada_legible(F)}")
    print("Veredicto:", _TEXTO[veredicto], "·", detalle)
    print("Recuerda escribir la respuesta con + C. Ángulos en radianes.")


widgets.interact(calculadora_verificador,
    f_txt=widgets.Text(value="cos(2*x)", description="f(x) ="),
    F_txt=widgets.Text(value="sin(2*x)", description="tu F(x) ="));

# %%
# Casos de prueba de la sección 2 (resultado conocido)
PRUEBAS_2 = [
    ("3x^2 con F = x^3 + 7: correcta", lambda: verificar(3 * x**2, x**3 + 7)[0] == "correcta"),
    ("cos 2x con F = sin 2x: incorrecta", lambda: verificar(sp.cos(2 * x), sp.sin(2 * x))[0] == "incorrecta"),
    ("1/x con F = ln|x| + 1: correcta", lambda: verificar(1 / x, sp.log(sp.Abs(x)) + 1)[0] == "correcta"),
    ("1/x con F = ln x: correcta solo en parte del dominio", lambda: verificar(1 / x, sp.log(x))[0] == "parcial"),
    ("(x+1)^2 con F = (x+1)^3/3: correcta y difiere de sympy en 1/3",
     lambda: verificar((x + 1)**2, (x + 1)**3 / 3)[0] == "correcta" and "1/3" in verificar((x + 1)**2, (x + 1)**3 / 3)[1]),
    ("1/cos^2 x con F = tan x: correcta", lambda: verificar(1 / sp.cos(x)**2, sp.tan(x))[0] == "correcta"),
]
for nombre, prueba in PRUEBAS_2:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 2).** Calcula a mano $\int\big(e^{2x}+\frac1x\big)dx$ y escribe tu antiderivada.

# %%
mi_F2 = None     # escribe tu antiderivada como texto, sin + C

if mi_F2 is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    f = sp.exp(2 * x) + 1 / x
    try:
        v, d = verificar(f, parsear(str(mi_F2)))
    except ValueError as err:
        v, d = "incorrecta", str(err)
    print("La calculadora da:", antiderivada(f), "+ C")
    print("Tu respuesta:", _TEXTO[v], "·", d)
    if v != "correcta":
        print("Revisa el 1/2 de la regla del argumento lineal y el valor absoluto del logaritmo.")
