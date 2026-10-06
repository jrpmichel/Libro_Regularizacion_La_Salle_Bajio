# ID: DIF-U1-NB00
# Notebook: dif/u1_funciones.ipynb · utilidades compartidas
# Repositorio: dif/u1_funciones/00_utilidades.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# # Dif. U1 · Funciones
# **Libro de regularización en matemáticas para ingeniería · ULSA Bajío**
#
# Este notebook acompaña la unidad 1 del bloque de cálculo diferencial. Tiene una sección por subtema, en el mismo orden que el libro:
#
# 1. Qué es una función: tablas de valores $(x,y)$
# 2. Función como regla de asignación y como gráfica
# 3. Dominio y rango
# 4. Álgebra de apoyo sobre la marcha
# 5. Funciones lineal y cuadrática
# 6. Funciones exponencial, logarítmica y trigonométricas
#
# **Cómo usarlo.** Ejecuta las celdas en orden (Entorno de ejecución → Ejecutar todas). Cada sección termina con una celda **Contrasta**: antes de ejecutarla, calcula a mano el caso que se indica y escribe tu resultado. La calculadora hace lo que le dices, no lo que querías decirle; por eso nunca la uses sin haber comprobado un caso que ya sabes resolver.
#
# **Etiqueta del repositorio.** `DIF-U1` (bloque, unidad). Las celdas de código de cada sección vienen de `dif/u1_funciones/NN_*.py`; los fragmentos que el libro imprime, de `dif/u1_funciones/fragmentos/`.
#
# **Convenciones.** Punto decimal. Los valores se *redondean* (no se truncan) a las cifras significativas que elijas. Para escribir expresiones usa `x` como variable, `*` para multiplicar, `**` o `^` para potencias y `sqrt(...)` para raíces.

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
x, y, h, a = sp.symbols("x y h a", real=True)

_LOCAL = {"x": x, "h": h, "a": a, "sqrt": sp.sqrt, "abs": sp.Abs, "exp": sp.exp,
          "log": sp.log, "pi": sp.pi, "sin": sp.sin, "cos": sp.cos, "tan": sp.tan}
_GLOBAL = {"Integer": sp.Integer, "Float": sp.Float, "Rational": sp.Rational,
           "Symbol": sp.Symbol, "Function": sp.Function}
_TRANSF = standard_transformations + (implicit_multiplication_application, convert_xor)


def parsear(texto, variables=("x",)):
    """Convierte un texto en expresión de sympy o explica qué salió mal."""
    if not isinstance(texto, str) or not texto.strip():
        raise ValueError("Escribe una expresión, por ejemplo: x**2 - 3*x + 1")
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
    libres = {str(s) for s in expr.free_symbols} - set(variables)
    if libres:
        raise ValueError(f"Símbolos no permitidos aquí: {sorted(libres)}. Solo puedes usar {list(variables)}.")
    return expr


def cifras(valor, n=4):
    """Redondea (no trunca) a n cifras significativas y devuelve texto."""
    v = float(valor)
    if not math.isfinite(v):
        return str(v)
    if v == 0:
        return "0"
    return np.format_float_positional(v, precision=n, unique=False, fractional=False, trim="-")


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


def comparar(mi_resultado, esperado):
    """Celda Contrasta: compara tu resultado a mano con el de la calculadora."""
    if mi_resultado is None or str(mi_resultado).strip() == "":
        print("Falta tu cálculo a mano. Escríbelo arriba y vuelve a ejecutar la celda.")
        return None
    try:
        mio = parsear(str(mi_resultado), variables=("x", "h", "a"))
        ok = sp.simplify(mio - sp.sympify(esperado)) == 0
    except ValueError as err:
        print(err); return None
    print("Coinciden." if ok else "No coinciden: revisa en qué paso te desviaste "
          "(¿sustituiste en todas las apariciones de x? ¿pusiste paréntesis?).")
    return ok
