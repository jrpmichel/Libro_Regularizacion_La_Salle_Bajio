# ID: DIF-U3-NB00
# Notebook: dif/u3_transformaciones.ipynb · utilidades compartidas
# Repositorio: dif/u3_transformaciones/00_utilidades.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# # Dif. U3 · Transformaciones gráficas
# **Libro de regularización en matemáticas para ingeniería · ULSA Bajío**
#
# Este notebook acompaña la unidad 3 del bloque de cálculo diferencial. Tiene una sección por subtema, en el mismo orden que el libro:
#
# 1. La forma general $y=a\,f\big(b(x-h)\big)+k$
# 2. Traslaciones y reflexiones
# 3. Estiramientos y compresiones
# 4. Composición de transformaciones
#
# **Cómo usarlo.** Ejecuta las celdas en orden (Entorno de ejecución → Ejecutar todas). Cada sección termina con una celda **Contrasta**: antes de ejecutarla, calcula a mano el caso que se indica y escribe tu resultado. Para saber hacia dónde se mueve una gráfica, sigue un punto: es más confiable que la intuición.
#
# **Etiqueta del repositorio.** `DIF-U3` (bloque, unidad). Las celdas de código de cada sección vienen de `dif/u3_transformaciones/NN_*.py`; los fragmentos que el libro imprime, de `dif/u3_transformaciones/fragmentos/`.
#
# **Convenciones.** Punto decimal. Los valores se *redondean* (no se truncan) a las cifras significativas que elijas. Para escribir expresiones usa `x` como variable, `*` para multiplicar, `**` o `^` para potencias, `sqrt(...)` para raíces, `abs(...)` para el valor absoluto, `exp(...)` para $e^x$ y `log(...)` para el logaritmo natural. Los ángulos de `sin`, `cos` y `tan` van en **radianes**.

# %% [markdown]
# ## 0. Versiones e imports

# %%
# Versiones mínimas probadas. En Colab ya vienen instaladas; esta celda solo lo asegura.
%pip install -q "numpy>=1.26" "sympy>=1.12" "matplotlib>=3.8" "ipywidgets>=8"

import sys, math, re
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sympy as sp
import ipywidgets as widgets
from IPython.display import display
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                        implicit_multiplication_application, convert_xor)
from sympy.calculus.util import continuous_domain, function_range

print("Python", sys.version.split()[0], "| numpy", np.__version__, "| sympy", sp.__version__,
      "| matplotlib", matplotlib.__version__, "| ipywidgets", widgets.__version__)

# %%
# Símbolos y utilidades compartidas por todo el notebook
x = sp.symbols("x", real=True)

_LOCAL = {"x": x, "sqrt": sp.sqrt, "abs": sp.Abs, "exp": sp.exp, "log": sp.log,
          "pi": sp.pi, "sin": sp.sin, "cos": sp.cos, "tan": sp.tan, "E": sp.E}
_GLOBAL = {"Integer": sp.Integer, "Float": sp.Float, "Rational": sp.Rational,
           "Symbol": sp.Symbol, "Function": sp.Function}
_TRANSF = standard_transformations + (implicit_multiplication_application, convert_xor)


def parsear(texto):
    """Convierte un texto en expresión de sympy en la variable x, o explica qué salió mal."""
    if not isinstance(texto, str) or not texto.strip():
        raise ValueError("Escribe una expresión, por ejemplo: (x**3 - 1)/(x - 1)")
    desconocidos = sorted(set(re.findall(r"[A-Za-z_]\w*", texto)) - set(_LOCAL))
    if desconocidos:
        raise ValueError(f"No reconozco {desconocidos}. Usa x como variable y, si hace falta, "
                         "sqrt, abs, exp, log (natural), sin, cos, tan o pi.")
    try:
        expr = parse_expr(texto, local_dict=_LOCAL, global_dict=dict(_GLOBAL),
                          transformations=_TRANSF)
    except Exception as err:
        raise ValueError("No pude leer la expresión. Usa x como variable, * para multiplicar, "
                         "** o ^ para potencias y sqrt(...) para raíces.") from err
    if expr.atoms(sp.core.function.AppliedUndef):
        raise ValueError("Hay una función que no conozco. Las permitidas son sqrt, abs, exp, log, sin, cos y tan.")
    libres = {str(s) for s in expr.free_symbols} - {"x"}
    if libres:
        raise ValueError(f"Símbolos no permitidos: {sorted(libres)}. Solo puedes usar x.")
    return expr


def numero(texto):
    """Lee un número escrito como 2, -0.5, 1/3, pi/2 o sqrt(2)."""
    try:
        v = sp.sympify(str(texto).replace("^", "**"), locals={"pi": sp.pi, "sqrt": sp.sqrt, "E": sp.E})
    except (sp.SympifyError, TypeError):
        raise ValueError(f"'{texto}' no es un número.") from None
    if v.free_symbols or not v.is_real:
        raise ValueError(f"'{texto}' no es un número real.")
    return v


def cifras(valor, n=4):
    """Redondea (no trunca) a n cifras significativas y devuelve texto."""
    if valor is None:
        return "no existe"
    if valor in (sp.oo, sp.S.Infinity) or (isinstance(valor, float) and valor == math.inf):
        return "+∞"
    if valor in (-sp.oo, sp.S.NegativeInfinity) or (isinstance(valor, float) and valor == -math.inf):
        return "-∞"
    v = float(valor)
    if math.isnan(v):
        return "no definida"
    if v == 0:
        return "0"
    return np.format_float_positional(v, precision=n, unique=False, fractional=False, trim="-")


def a_texto(valor):
    """Límite de sympy en texto: número exacto si es racional, ±∞ o 'no existe'."""
    if valor is None or valor is sp.nan:
        return "no existe"
    if valor == sp.oo:
        return "+∞"
    if valor == -sp.oo:
        return "-∞"
    if valor == sp.zoo:
        return "∞ sin signo (no existe)"
    if isinstance(valor, sp.AccumBounds):
        return "no existe (oscila)"
    return str(valor) if valor.is_Rational else f"{valor} ≈ {cifras(valor, 6)}"


def conjunto_a_texto(s):
    """Intervalos de sympy en notación de intervalos."""
    if s == sp.S.Reals:
        return "(-∞, ∞) = ℝ"
    if s == sp.S.EmptySet:
        return "∅ (vacío)"
    if isinstance(s, sp.Union):
        return " ∪ ".join(conjunto_a_texto(t) for t in sorted(s.args, key=lambda t: float(t.inf)))
    if isinstance(s, sp.Interval):
        f = lambda p: ("∞" if p == sp.oo else "-∞" if p == -sp.oo
                       else str(p) if p.is_Rational else "≈" + cifras(float(p), 4))
        return ("(" if s.left_open else "[") + f(s.inf) + ", " + f(s.sup) + (")" if s.right_open else "]")
    return str(s)


def evaluar(f, v):
    """f(v) en punto flotante; nan si f no está definida en v."""
    with np.errstate(all="ignore"):
        try:
            r = complex(f(v))
        except (ZeroDivisionError, ValueError, OverflowError):
            return math.nan
    return r.real if abs(r.imag) < 1e-12 else math.nan


def cerca(u, v, tol=1e-9):
    return abs(u - v) < tol


def _rechaza(fn):
    """True si fn() lanza ValueError (para probar que una calculadora rechaza una entrada)."""
    try:
        fn()
    except ValueError:
        return True
    return False


# Funciones base del libro, con tres puntos clave de cada una (x0, y0)
BASES = {
    "x²":    (x**2,        [(-1, 1), (0, 0), (2, 4)]),
    "|x|":   (sp.Abs(x),   [(-1, 1), (0, 0), (2, 2)]),
    "√x":    (sp.sqrt(x),  [(0, 0), (1, 1), (4, 2)]),
    "1/x":   (1/x,         [(-1, -1), (1, 1), (2, sp.Rational(1, 2))]),
    "e^x":   (sp.exp(x),   [(-1, sp.exp(-1)), (0, 1), (1, sp.E)]),
    "sen x": (sp.sin(x),   [(0, 0), (sp.pi/2, 1), (sp.pi, 0)]),
}


def transformar(f, a, b, h, k):
    """g(x) = a f(b(x - h)) + k. Rechaza a = 0 o b = 0, que aplastan la gráfica."""
    if a == 0:
        raise ValueError("a = 0 convierte la gráfica en la recta y = k: ya no es una transformación de f.")
    if b == 0:
        raise ValueError("b = 0 convierte la función en la constante a·f(0) + k (o en nada, si f(0) no existe): "
                         "ya no es una transformación de f.")
    a, b, h, k = (sp.nsimplify(v) for v in (a, b, h, k))
    return a * f.subs(x, b * (x - h)) + k


def imagen(punto, a, b, h, k):
    """(x0, y0) de la gráfica de f  ->  (h + x0/b, k + a·y0) de la gráfica de g."""
    x0, y0 = punto
    return (sp.nsimplify(h) + sp.sympify(x0) / sp.nsimplify(b), sp.nsimplify(k) + sp.nsimplify(a) * sp.sympify(y0))
