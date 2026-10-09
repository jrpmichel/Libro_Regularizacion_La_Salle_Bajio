# ID: DIF-U4-NB00
# Notebook: dif/u4_derivada.ipynb · utilidades compartidas
# Repositorio: dif/u4_derivada/00_utilidades.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# # Dif. U4 · La derivada
# **Libro de regularización en matemáticas para ingeniería · ULSA Bajío**
#
# Este notebook acompaña la unidad 4 del bloque de cálculo diferencial. Tiene una sección por subtema, en el mismo orden que el libro:
#
# 1. Pendiente de una recta y pendiente de la secante
# 2. Paso al límite: definición de derivada
# 3. Primeros ejemplos: $x^2$, $x^3$, $1/x$, $\sqrt{x}$
# 4. Notaciones $f'(x)$ y $dy/dx$; recta tangente
# 5. La derivada como razón de cambio
# 6. Propiedades: derivable implica continua; puntos no derivables
#
# **Cómo usarlo.** Ejecuta las celdas en orden (Entorno de ejecución → Ejecutar todas). Cada sección termina con una celda **Contrasta**: antes de ejecutarla, calcula a mano el caso que se indica y escribe tu resultado. En la definición de derivada, $h$ es una variable que tiende a cero, nunca un número fijo: si tu resultado todavía depende de $h$, falta el paso al límite.
#
# **Etiqueta del repositorio.** `DIF-U4` (bloque, unidad). Las celdas de código de cada sección vienen de `dif/u4_derivada/NN_*.py`; los fragmentos que el libro imprime, de `dif/u4_derivada/fragmentos/`.
#
# **Convenciones.** Punto decimal. Los valores se *redondean* (no se truncan) a las cifras significativas que elijas. Para escribir expresiones usa `x` como variable (también cuando represente el tiempo), `*` para multiplicar, `**` o `^` para potencias, `sqrt(...)` para raíces cuadradas, `cbrt(...)` para la raíz cúbica real, `abs(...)` para el valor absoluto, `exp(...)` para $e^x$ y `log(...)` para el logaritmo natural. Los ángulos de `sin`, `cos` y `tan` van en **radianes**.

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
                                        implicit_multiplication_application, convert_xor, rationalize)
from sympy.calculus.util import continuous_domain, function_range

print("Python", sys.version.split()[0], "| numpy", np.__version__, "| sympy", sp.__version__,
      "| matplotlib", matplotlib.__version__, "| ipywidgets", widgets.__version__)

# %%
# Símbolos y utilidades compartidas por todo el notebook
x = sp.symbols("x", real=True)
h = sp.symbols("h", real=True)     # incremento: variable que tiende a cero

_LOCAL = {"x": x, "sqrt": sp.sqrt, "abs": sp.Abs, "exp": sp.exp, "log": sp.log,
          "cbrt": lambda e: sp.sign(e) * sp.Abs(e)**sp.Rational(1, 3),
          "pi": sp.pi, "sin": sp.sin, "cos": sp.cos, "tan": sp.tan, "E": sp.E}
_GLOBAL = {"Integer": sp.Integer, "Float": sp.Float, "Rational": sp.Rational,
           "Symbol": sp.Symbol, "Function": sp.Function}
# rationalize: 0.1 se lee como 1/10, para que los límites salgan exactos
_TRANSF = standard_transformations + (implicit_multiplication_application, convert_xor, rationalize)


def parsear(texto):
    """Convierte un texto en expresión de sympy en la variable x, o explica qué salió mal."""
    if not isinstance(texto, str) or not texto.strip():
        raise ValueError("Escribe una expresión, por ejemplo: (x**3 - 1)/(x - 1)")
    desconocidos = sorted(set(re.findall(r"[A-Za-z_]\w*", texto)) - set(_LOCAL))
    if desconocidos:
        raise ValueError(f"No reconozco {desconocidos}. Usa x como variable y, si hace falta, "
                         "sqrt, cbrt, abs, exp, log (natural), sin, cos, tan o pi.")
    try:
        expr = parse_expr(texto, local_dict=_LOCAL, global_dict=dict(_GLOBAL),
                          transformations=_TRANSF)
    except Exception as err:
        raise ValueError("No pude leer la expresión. Usa x como variable, * para multiplicar, "
                         "** o ^ para potencias y sqrt(...) para raíces.") from err
    if expr.atoms(sp.core.function.AppliedUndef):
        raise ValueError("Hay una función que no conozco. Las permitidas son sqrt, cbrt, abs, exp, log, sin, cos y tan.")
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
    if valor.has(sp.I):
        return "no definido (fuera del dominio)"
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


def definida_en(f, a):
    """True si la expresión f tiene valor real finito en x = a."""
    v = f.subs(x, a)
    return bool(v.is_real) and bool(v.is_finite)


def pendiente_secante(f, a, b):
    """(f(b) - f(a)) / (b - a), exacta. Rechaza a = b y puntos fuera del dominio."""
    a, b = sp.nsimplify(a), sp.nsimplify(b)
    if a == b:
        raise ValueError("a y b son iguales: con un solo punto no hay secante (Δx = 0).")
    for p in (a, b):
        if not definida_en(f, p):
            raise ValueError(f"f no está definida en x = {p}: ese punto no está en la gráfica.")
    return sp.simplify((f.subs(x, b) - f.subs(x, a)) / (b - a))


def cociente(f, a):
    """Cociente incremental (f(a + h) - f(a)) / h como expresión en h."""
    a = sp.nsimplify(a)
    if not definida_en(f, a):
        raise ValueError(f"f no está definida en x = {a}: no hay punto de tangencia.")
    return (f.subs(x, a + h) - f.subs(x, a)) / h


def lado_definido(f, a, signo):
    """True si f tiene valores reales justo a la izquierda (signo = -1) o a la derecha (signo = +1) de a."""
    v = f.subs(x, sp.nsimplify(a) + signo * sp.Rational(1, 10**6))
    return bool(v.is_real) and bool(v.is_finite)


def lado_texto(valor):
    """Derivada lateral en texto; None significa que ese lado queda fuera del dominio."""
    return "no definida (ese lado está fuera del dominio)" if valor is None else a_texto(valor)


def derivada_en(f, a):
    """(izquierda, derecha) del límite del cociente cuando h -> 0, y f'(a) o None si no existe.
    Un lateral vale None si f no está definida de ese lado de a (extremo del dominio, como √x en 0)."""
    q = cociente(f, a)
    izq = sp.limit(q, h, 0, "-") if lado_definido(f, a, -1) else None
    der = sp.limit(q, h, 0, "+") if lado_definido(f, a, +1) else None
    existe = (izq is not None and der is not None and izq.is_finite and der.is_finite
              and sp.simplify(izq - der) == 0)
    return izq, der, (sp.simplify(der) if existe else None)
