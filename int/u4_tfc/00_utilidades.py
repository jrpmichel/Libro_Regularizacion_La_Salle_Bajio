# ID: INT-U4-NB00
# Notebook: int/u4_tfc.ipynb · utilidades compartidas
# Repositorio: int/u4_tfc/00_utilidades.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# # Int. U4 · Teorema fundamental del cálculo
# **Libro de regularización en matemáticas para ingeniería · ULSA Bajío**
#
# Este notebook acompaña la unidad 4 del bloque de cálculo integral. Tiene una sección por subtema, en el mismo orden que el libro:
#
# 1. Primera parte: la derivada de la integral
# 2. Segunda parte: evaluación con la antiderivada
# 3. Evaluación de integrales definidas (integral con signo y área total)
# 4. Teorema del valor medio para integrales
# 5. Comprobación numérica contra la suma de Riemann
#
# **Cómo usarlo.** Ejecuta las celdas en orden (Entorno de ejecución → Ejecutar todas). Cada sección termina con una celda **Contrasta**: antes de ejecutarla, resuelve a mano el caso que se indica y escribe tu resultado. La referencia de la calculadora aparece solo después de que escribes el tuyo.
#
# **Etiqueta del repositorio.** `INT-U4` (bloque, unidad). Las celdas de código de cada sección vienen de `int/u4_tfc/NN_*.py`; los fragmentos que el libro imprime, de `int/u4_tfc/fragmentos/`.
#
# **Convenciones.** Punto decimal. Los valores se *redondean* (no se truncan) a las cifras significativas que elijas. Para escribir funciones usa `x` como variable, `*` para multiplicar, `**` o `^` para potencias, `sqrt(...)`, `exp(...)`, `log(...)` (natural) y `abs(...)`. Los ángulos de `sin`, `cos`, `tan` y `sec` van en **radianes**. Antes de aplicar $F(b)-F(a)$, las calculadoras revisan que $f$ sea continua en $[a,b]$.

# %% [markdown]
# ## 0. Versiones e imports

# %%
# Versiones mínimas probadas. En Colab ya vienen instaladas; esta celda solo lo asegura.
%pip install -q "numpy>=1.26" "sympy>=1.12" "matplotlib>=3.8" "ipywidgets>=8" "scipy>=1.11"

import sys, math, re
from decimal import Decimal, ROUND_HALF_UP
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sympy as sp
import ipywidgets as widgets
from IPython.display import display
from scipy.optimize import brentq
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                        implicit_multiplication_application, convert_xor, rationalize)

import scipy
print("Python", sys.version.split()[0], "| numpy", np.__version__, "| scipy", scipy.__version__, "| sympy", sp.__version__,
      "| matplotlib", matplotlib.__version__, "| ipywidgets", widgets.__version__)

# %%
# Símbolos y utilidades compartidas por todo el notebook
x = sp.symbols("x", real=True)

_FUNCIONES = {"sqrt": sp.sqrt, "abs": sp.Abs, "exp": sp.exp, "log": sp.log, "ln": sp.log,
              "pi": sp.pi, "sin": sp.sin, "sen": sp.sin, "cos": sp.cos, "tan": sp.tan, "sec": sp.sec,
              "E": sp.E, "e": sp.E}
_GLOBAL = {"Integer": sp.Integer, "Float": sp.Float, "Rational": sp.Rational,
           "Symbol": sp.Symbol, "Function": sp.Function}
# rationalize: 0.1 se lee como 1/10, para que los resultados salgan exactos
_TRANSF = standard_transformations + (implicit_multiplication_application, convert_xor, rationalize)
_OTRAS_VARIABLES = ("t", "u", "w", "q", "s")      # se leen como x si son la única variable
_C = sp.Symbol("C")


def parsear(texto, variables=("x",)):
    """Convierte un texto en una expresión de sympy o explica qué salió mal.
    Acepta sen y ln, la letra e como número de Euler, una constante + C (se descarta)
    y otra letra (t, u, w, q o s) en lugar de x si es la única variable."""
    if not isinstance(texto, str) or not texto.strip():
        raise ValueError("Escribe una función, por ejemplo: x**2")
    nombres = set(re.findall(r"[A-Za-z_]\w*", texto))
    otras = [v for v in _OTRAS_VARIABLES if v in nombres]
    if "x" not in nombres and len(otras) == 1:
        texto = re.sub(rf"(?<![A-Za-z_]){otras[0]}(?![A-Za-z_0-9])", "x", texto)
        nombres = set(re.findall(r"[A-Za-z_]\w*", texto))
    local = {**_FUNCIONES, "x": x, "C": _C}
    desconocidos = sorted(nombres - set(local))
    if desconocidos:
        raise ValueError(f"No reconozco {desconocidos}. Usa x (o t) como variable y, si hace falta, "
                         "sqrt, exp, log o ln, abs, sin o sen, cos, tan o pi.")
    try:
        expr = parse_expr(texto, local_dict=local, global_dict=dict(_GLOBAL), transformations=_TRANSF)
    except Exception as err:
        raise ValueError("No pude leer la función. Usa * para multiplicar, ** o ^ para potencias "
                         "y sqrt(...) para raíces.") from err
    if not isinstance(expr, sp.Expr) or expr.atoms(sp.core.function.AppliedUndef):
        raise ValueError("Escribe una expresión sin el signo = y solo con las funciones permitidas.")
    expr = expr.subs(_C, 0)                       # la constante de integración no cambia la derivada
    if expr.free_symbols - {x}:
        raise ValueError("Usa una sola variable.")
    return expr


def lista_numeros(texto, nombre="la lista"):
    """'0, 1.5, 3' -> [0.0, 1.5, 3.0], o un mensaje claro si algo no es número."""
    partes = [p for p in re.split(r"[,\s;]+", str(texto).strip()) if p]
    if not partes:
        raise ValueError(f"Escribe {nombre} como números separados por comas, por ejemplo: 0, 1, 3, 4")
    try:
        return [float(p) for p in partes]
    except ValueError as err:
        raise ValueError(f"En {nombre} hay algo que no es un número: usa punto decimal y comas entre valores.") from err


def cifras(valor, n=4):
    """Redondea (no trunca) a n cifras significativas con la regla escolar (5 sube) y conserva los ceros finales."""
    v = float(valor)
    if math.isnan(v):
        return "no definido"
    if math.isinf(v):
        return "+∞" if v > 0 else "-∞"
    if v == 0:
        return "0"
    d = Decimal(repr(v))
    q = d.quantize(Decimal(1).scaleb(d.adjusted() - n + 1), rounding=ROUND_HALF_UP)
    if q.adjusted() != d.adjusted():          # 9.996 con 3 cifras pasa a 10.0: se reajusta la posición
        q = d.quantize(Decimal(1).scaleb(q.adjusted() - n + 1), rounding=ROUND_HALF_UP)
    if q.adjusted() >= n:                     # enteros grandes: sin notación científica
        return f"{int(q):d}"
    return format(q, "f")


def cerca(p, q, tol=1e-9):
    return abs(float(p) - float(q)) < tol


def _rechaza(fn):
    """True si fn() lanza ValueError (para probar que una calculadora rechaza una entrada)."""
    try:
        fn()
    except ValueError:
        return True
    return False


def valor_real(f, c):
    """f(c) como float, o ValueError si no es un número real finito."""
    v = complex(sp.N(f.subs(x, c)))
    if abs(v.imag) > 1e-12 or not math.isfinite(v.real):
        raise ValueError(f"f no tiene valor real en x = {cifras(float(c), 6)}: revisa el dominio.")
    return v.real


def numerica(f):
    """Versión numpy de f que siempre devuelve un arreglo del tamaño de la entrada."""
    g = sp.lambdify(x, f, "numpy")
    return lambda t: np.asarray(g(np.asarray(t, dtype=float)), dtype=float) * np.ones_like(np.asarray(t, dtype=float))


# Antiderivadas (las mismas herramientas de la Unidad 3)
_PERMITIDAS = (sp.exp, sp.log, sp.sin, sp.cos, sp.tan, sp.sec, sp.Abs, sp.atan, sp.asin, sp.acos)
_PUNTOS = [-2.9, -2.1, -1.3, -0.7, 0.7, 1.3, 2.1, 2.9]


def con_valor_absoluto(F):
    """ln(u) -> ln|u| salvo que u sea siempre positivo: la antiderivada de u'/u vale donde u > 0 y donde u < 0."""
    return F.replace(sp.log, lambda a: sp.log(a) if (a.is_positive or a.has(sp.Abs)) else sp.log(sp.Abs(a)))


def sin_abs_en_log(f):
    """ln|u| -> ln u, solo para integrar (la derivada de las dos es u'/u); el valor absoluto se repone al final."""
    return f.replace(lambda e: isinstance(e, sp.log) and e.args[0].has(sp.Abs),
                     lambda e: sp.log(e.args[0].replace(sp.Abs, lambda a: a)))


def misma_funcion(g, h):
    """True si g y h coinciden: primero con sympy; si no decide, en los puntos donde las dos son reales."""
    d = sp.simplify(g - h)
    if d == 0:
        return True
    comparados = 0
    for c in _PUNTOS:
        try:
            if abs(valor_real(g, c) - valor_real(h, c)) > 1e-9:
                return False
            comparados += 1
        except (ValueError, TypeError, ZeroDivisionError):
            continue
    return comparados >= 3


def antiderivada(f):
    """Una antiderivada de f (con C = 0), simplificada y con ln|..| donde corresponde."""
    F = sp.integrate(sin_abs_en_log(f), x)
    if F.has(sp.Integral):
        raise ValueError("sympy no encontró una antiderivada para esta función.")
    raras = {type(a).__name__ for a in F.atoms(sp.Function) if not isinstance(a, _PERMITIDAS)}
    if raras:
        raise ValueError(f"la antiderivada usa funciones especiales ({', '.join(sorted(raras))}): no es elemental. "
                         "Calcula la integral definida con las sumas de la sección 5 de este notebook.")
    F = sp.expand(F) if F.is_polynomial(x) else con_valor_absoluto(sp.trigsimp(sp.simplify(F)))
    if not misma_funcion(derivada_legible(F), sin_abs_en_log(f)):
        raise ValueError("la comprobación F' = f falló; revisa la función.")
    return F


def derivada_legible(F):
    """F' sin el Piecewise que sympy agrega al derivar ln|u| (la derivada de ln|u| es u'/u)."""
    G = F.replace(lambda e: isinstance(e, sp.log) and e.args[0].has(sp.Abs),
                  lambda e: sp.log(e.args[0].replace(sp.Abs, lambda a: a)))
    return sp.simplify(sp.diff(G, x))


def continua_en(f, a, b):
    """True si f es continua en todo el intervalo cerrado entre a y b."""
    lo, hi = sorted((sp.nsimplify(a), sp.nsimplify(b)))
    dom = sp.calculus.util.continuous_domain(f, x, sp.Interval(lo, hi))
    return dom == sp.Interval(lo, hi)
