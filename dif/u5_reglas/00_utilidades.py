# ID: DIF-U5-NB00
# Notebook: dif/u5_reglas.ipynb · utilidades compartidas
# Repositorio: dif/u5_reglas/00_utilidades.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# # Dif. U5 · Reglas de derivación
# **Libro de regularización en matemáticas para ingeniería · ULSA Bajío**
#
# Este notebook acompaña la unidad 5 del bloque de cálculo diferencial. Tiene una sección por subtema, en el mismo orden que el libro:
#
# 1. Reglas de constante, potencia y suma
# 2. Reglas de producto y cociente
# 3. Regla de la cadena
# 4. Derivadas de orden superior
# 5. Derivadas de exponencial, logaritmo y funciones trigonométricas
# 6. Derivación implícita
# 7. Derivación logarítmica
#
# **Cómo usarlo.** Ejecuta las celdas en orden (Entorno de ejecución → Ejecutar todas). Cada sección termina con una celda **Contrasta**: antes de ejecutarla, deriva a mano el caso que se indica y escribe tu resultado. Comprueba siempre una derivada en al menos dos puntos y evita los que simplifican demasiado ($x=0$, $x=1$): en ellos una fórmula equivocada puede coincidir con la correcta.
#
# **Etiqueta del repositorio.** `DIF-U5` (bloque, unidad). Las celdas de código de cada sección vienen de `dif/u5_reglas/NN_*.py`; los fragmentos que el libro imprime, de `dif/u5_reglas/fragmentos/`.
#
# **Convenciones.** Punto decimal. Los valores se *redondean* (no se truncan) a las cifras significativas que elijas. Para escribir expresiones usa `x` como variable, `*` para multiplicar, `**` o `^` para potencias, `sqrt(...)` para raíces cuadradas, `exp(...)` para $e^x$ y `log(...)` para el logaritmo natural. Los ángulos de `sin`, `cos` y `tan` van en **radianes**, salvo en la calculadora de la sección 5, que tiene un selector.

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

print("Python", sys.version.split()[0], "| numpy", np.__version__, "| sympy", sp.__version__,
      "| matplotlib", matplotlib.__version__, "| ipywidgets", widgets.__version__)

# %%
# Símbolos y utilidades compartidas por todo el notebook
x, y, u = sp.symbols("x y u", real=True)

_FUNCIONES = {"sqrt": sp.sqrt, "abs": sp.Abs, "exp": sp.exp, "log": sp.log,
              "pi": sp.pi, "sin": sp.sin, "cos": sp.cos, "tan": sp.tan, "E": sp.E}
_GLOBAL = {"Integer": sp.Integer, "Float": sp.Float, "Rational": sp.Rational,
           "Symbol": sp.Symbol, "Function": sp.Function}
# rationalize: 0.1 se lee como 1/10, para que las derivadas salgan exactas
_TRANSF = standard_transformations + (implicit_multiplication_application, convert_xor, rationalize)
_SIMBOLOS = {"x": x, "y": y, "u": u}


def parsear(texto, variables=("x",), funciones=None):
    """Convierte un texto en expresión de sympy o explica qué salió mal.
    variables: nombres permitidos (x, y o u). funciones: reemplaza el diccionario de funciones."""
    if not isinstance(texto, str) or not texto.strip():
        raise ValueError("Escribe una expresión, por ejemplo: x**3 - 3*x")
    fun = dict(_FUNCIONES if funciones is None else funciones)
    local = {**fun, **{v: _SIMBOLOS[v] for v in variables}}
    desconocidos = sorted(set(re.findall(r"[A-Za-z_]\w*", texto)) - set(local))
    if desconocidos:
        raise ValueError(f"No reconozco {desconocidos}. Usa {', '.join(variables)} como variable y, si hace falta, "
                         "sqrt, exp, log (natural), sin, cos, tan o pi.")
    try:
        expr = parse_expr(texto, local_dict=local, global_dict=dict(_GLOBAL), transformations=_TRANSF)
    except Exception as err:
        raise ValueError("No pude leer la expresión. Usa * para multiplicar, ** o ^ para potencias "
                         "y sqrt(...) para raíces.") from err
    if not isinstance(expr, sp.Expr):
        raise ValueError("Escribe una expresión, no una ecuación (sin el signo =).")
    if expr.atoms(sp.core.function.AppliedUndef):
        raise ValueError("Hay una función que no conozco. Las permitidas son sqrt, exp, log, sin, cos y tan.")
    libres = {str(s) for s in expr.free_symbols} - set(variables)
    if libres:
        raise ValueError(f"Símbolos no permitidos: {sorted(libres)}. Solo puedes usar {list(variables)}.")
    return expr


def cifras(valor, n=4):
    """Redondea (no trunca) a n cifras significativas y devuelve texto."""
    v = float(valor)
    if math.isnan(v):
        return "no definido"
    if math.isinf(v):
        return "+∞" if v > 0 else "-∞"
    if v == 0:
        return "0"
    return np.format_float_positional(v, precision=n, unique=False, fractional=False, trim="-")


def iguales(e1, e2):
    """True si las dos expresiones son equivalentes."""
    return sp.simplify(sp.sympify(e1) - sp.sympify(e2)) == 0


def evaluar_en(expr, var, a):
    """Valor exacto de expr en var = a; ValueError si no es un número real finito."""
    v = sp.simplify(expr.subs(var, sp.nsimplify(a)))
    if not (v.is_real and v.is_finite):
        raise ValueError(f"la expresión no tiene valor real en {var} = {a} (fuera del dominio).")
    return v


def cociente_num(expr, var, a, hh=1e-6):
    """Cociente incremental numérico (f(a + h) - f(a))/h, como comprobación independiente."""
    f = sp.lambdify(var, expr, "math")
    return (f(float(a) + hh) - f(float(a))) / hh


def cerca(p, q, tol=1e-9):
    return abs(float(p) - float(q)) < tol


def _rechaza(fn):
    """True si fn() lanza ValueError (para probar que una calculadora rechaza una entrada)."""
    try:
        fn()
    except ValueError:
        return True
    return False


def contrasta(mi_texto, referencia, variables=("x",), pista=""):
    """Celda Contrasta: compara tu resultado a mano (texto) con el de la calculadora."""
    if mi_texto is None or str(mi_texto).strip() == "":
        print("Falta tu cálculo a mano. Escríbelo arriba y vuelve a ejecutar la celda.")
        return None
    try:
        mio = parsear(str(mi_texto), variables)
    except ValueError as err:
        print(err); return None
    ok = iguales(mio, referencia)
    print("Coinciden." if ok else "No coinciden. " + pista)
    return ok
