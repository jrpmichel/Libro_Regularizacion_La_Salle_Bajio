# ID: INT-U2-NB00
# Notebook: int/u2_riemann.ipynb · utilidades compartidas
# Repositorio: int/u2_riemann/00_utilidades.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# # Int. U2 · Sumas de Riemann
# **Libro de regularización en matemáticas para ingeniería · ULSA Bajío**
#
# Este notebook acompaña la unidad 2 del bloque de cálculo integral. Tiene una sección por subtema, en el mismo orden que el libro:
#
# 1. Suma izquierda, derecha y de punto medio
# 2. Notación sigma y fórmulas cerradas
# 3. La integral definida como límite de sumas
# 4. Propiedades básicas de la integral definida
# 5. Regla del trapecio y regla de Simpson, con tabla de error contra n
#
# **Cómo usarlo.** Ejecuta las celdas en orden (Entorno de ejecución → Ejecutar todas). Cada sección termina con una celda **Contrasta**: antes de ejecutarla, resuelve a mano el caso que se indica y escribe tu resultado. La referencia de la calculadora aparece solo después de que escribes el tuyo.
#
# **Etiqueta del repositorio.** `INT-U2` (bloque, unidad). Las celdas de código de cada sección vienen de `int/u2_riemann/NN_*.py`; los fragmentos que el libro imprime, de `int/u2_riemann/fragmentos/`.
#
# **Convenciones.** Punto decimal. Los valores se *redondean* (no se truncan) a las cifras significativas que elijas. Para escribir funciones usa `x` como variable (en las sumas, `k` es el índice y `n` el límite superior), `*` para multiplicar, `**` o `^` para potencias, `sqrt(...)`, `exp(...)` y `log(...)` (natural). Los ángulos de `sin`, `cos` y `tan` van en **radianes**.

# %% [markdown]
# ## 0. Versiones e imports

# %%
# Versiones mínimas probadas. En Colab ya vienen instaladas; esta celda solo lo asegura.
%pip install -q "numpy>=1.26" "sympy>=1.12" "matplotlib>=3.8" "ipywidgets>=8"

import sys, math, re
from decimal import Decimal, ROUND_HALF_UP
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sympy as sp
import ipywidgets as widgets
from IPython.display import display
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                        implicit_multiplication_application, convert_xor, rationalize)

print("Python", sys.version.split()[0], "| numpy", np.__version__, "| sympy", sp.__version__,
      "| matplotlib", matplotlib.__version__, "| ipywidgets", widgets.__version__)

# %%
# Símbolos y utilidades compartidas por todo el notebook
x = sp.symbols("x", real=True)
k, n = sp.symbols("k n", integer=True, positive=True)

_FUNCIONES = {"sqrt": sp.sqrt, "abs": sp.Abs, "exp": sp.exp, "log": sp.log,
              "pi": sp.pi, "sin": sp.sin, "cos": sp.cos, "tan": sp.tan, "E": sp.E}
_GLOBAL = {"Integer": sp.Integer, "Float": sp.Float, "Rational": sp.Rational,
           "Symbol": sp.Symbol, "Function": sp.Function}
# rationalize: 0.1 se lee como 1/10, para que los resultados salgan exactos
_TRANSF = standard_transformations + (implicit_multiplication_application, convert_xor, rationalize)


def parsear(texto, variables=("x",)):
    """Convierte un texto en una expresión de sympy o explica qué salió mal.
    variables: nombres permitidos (x, k o n)."""
    if not isinstance(texto, str) or not texto.strip():
        raise ValueError("Escribe una función, por ejemplo: x**2")
    simbolos = {"x": x, "k": k, "n": n}
    local = {**_FUNCIONES, **{v: simbolos[v] for v in variables}}
    desconocidos = sorted(set(re.findall(r"[A-Za-z_]\w*", texto)) - set(local))
    if desconocidos:
        raise ValueError(f"No reconozco {desconocidos}. Usa {', '.join(variables)} como variable y, si hace falta, "
                         "sqrt, exp, log (natural), sin, cos, tan o pi.")
    try:
        expr = parse_expr(texto, local_dict=local, global_dict=dict(_GLOBAL), transformations=_TRANSF)
    except Exception as err:
        raise ValueError("No pude leer la función. Usa * para multiplicar, ** o ^ para potencias "
                         "y sqrt(...) para raíces.") from err
    if not isinstance(expr, sp.Expr) or expr.atoms(sp.core.function.AppliedUndef):
        raise ValueError("Escribe una expresión sin el signo = y solo con las funciones permitidas.")
    if {str(v) for v in expr.free_symbols} - set(variables):
        raise ValueError(f"Solo puedes usar {', '.join(variables)}.")
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
