# ID: ALG-U1-NB00
# Notebook: alg/u1_vectores.ipynb · utilidades compartidas
# Repositorio: alg/u1_vectores/00_utilidades.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# # Álg. U1 · Vectores, rectas y planos
# **Libro de regularización en matemáticas para ingeniería · ULSA Bajío**
#
# Este notebook acompaña la unidad 1 del bloque de álgebra lineal. Tiene una sección por subtema, en el mismo orden que el libro:
#
# 1. Vectores en 2D y 3D: componentes y magnitud
# 2. Suma, resta y vector unitario
# 3. Producto punto y proyección
# 4. Producto cruz
# 5. Rectas y planos en el espacio
#
# **Cómo usarlo.** Ejecuta las celdas en orden (Entorno de ejecución → Ejecutar todas). Cada sección tiene una calculadora, sus casos de prueba y termina con una celda **Contrasta**: antes de ejecutarla, resuelve a mano el caso que se indica y escribe tu resultado. La referencia de la calculadora aparece solo después de que escribes el tuyo.
#
# **Etiqueta del repositorio.** `ALG-U1` (bloque, unidad). Las celdas de código de cada sección vienen de `alg/u1_vectores/NN_*.py`; los fragmentos que el libro imprime, de `alg/u1_vectores/fragmentos/`.
#
# **Convenciones.** Punto decimal. Un vector se escribe como números separados por comas, por ejemplo `1, -2, 3`; si quieres, enciérralo entre paréntesis o corchetes: `(1, -2, 3)` o `[1 -2 3]` (sin comas, separa con espacios). Cada componente acepta fracciones como `2/3`, raíces como `sqrt(2)` y `pi`. Los ángulos salen en **grados** o **radianes** según el control `unidad`, y la unidad siempre aparece junto al valor. Los valores se *redondean* (no se truncan) a las cifras significativas que elijas; los muy grandes o muy pequeños se escriben en notación científica, por ejemplo `1.234×10^-7`.

# %% [markdown]
# ## 0. Versiones e imports

# %%
# Versiones mínimas probadas. En Colab ya vienen instaladas; esta celda solo lo asegura.
%pip install -q "numpy>=1.26" "sympy>=1.12" "matplotlib>=3.8" "ipywidgets>=8"

import sys, math, re
from decimal import Decimal, ROUND_HALF_UP
from fractions import Fraction
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registra la proyección "3d")
import sympy as sp
from sympy.core.evalf import PrecisionExhausted
import ipywidgets as widgets
from IPython.display import display
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                        implicit_multiplication_application, convert_xor, rationalize)

print("Python", sys.version.split()[0], "| numpy", np.__version__, "| sympy", sp.__version__,
      "| matplotlib", matplotlib.__version__, "| ipywidgets", widgets.__version__)

# %%
# Utilidades compartidas por todo el notebook
_PERMITIDOS = {"sqrt": sp.sqrt, "pi": sp.pi}
_GLOBAL = {"Integer": sp.Integer, "Float": sp.Float, "Rational": sp.Rational,
           "Add": sp.Add, "Mul": sp.Mul, "Pow": sp.Pow}
# rationalize: 0.1 se lee como 1/10, para que las fracciones salgan exactas antes de pasar a decimal
_TRANSF = standard_transformations + (implicit_multiplication_application, convert_xor, rationalize)
_CARACTERES = re.compile(r"^[0-9A-Za-z.+\-*/^() ]+$")
_LITERAL = re.compile(r"\d*\.\d+(?:[eE][+-]?\d+)?|\d+\.?(?:[eE][+-]?\d+)?")
_MAXIMO, _MINIMO = 1e100, 1e-100     # fuera de este rango los productos desbordan el punto flotante
_CIERRE = {"(": ")", "[": "]", "<": ">", "⟨": "⟩"}


def _contrae(preposicion, nombre):
    """'de' + 'el vector' -> 'del vector';  'a' + 'el valor' -> 'al valor'."""
    return f"{preposicion}l {nombre[3:]}" if nombre.startswith("el ") else f"{preposicion} {nombre}"


def _revisa_rango(v, nombre):
    """Rechaza componentes tan grandes o tan pequeñas que los productos se desbordan."""
    if v != 0 and not (_MINIMO * (1 - 1e-9) <= abs(v) <= _MAXIMO * (1 + 1e-9)):
        raise ValueError(f"{nombre} tiene un número fuera del rango 1e-100 a 1e100 (en valor absoluto); "
                         "cambia de unidades, por ejemplo de mm a m.")


def _potencia_excesiva(p):
    """True si p es una potencia con exponente mayor que 1000 o con valor mayor que 10^2000:
    evaluarla de forma exacta tardaría demasiado (9^9^9)."""
    if not isinstance(p, sp.Pow):
        return False
    if not abs(sp.N(p.exp, 15)) <= 1000:
        return True
    v = sp.N(p, 15)
    return bool(v.is_finite) and abs(v) > sp.Integer(10)**2000


def _lee_numero(texto, nombre):
    """'2/3', 'sqrt(2)', '-1.5e3', '3pi' -> float, o un mensaje claro si algo no es un número real."""
    t = texto.strip().replace("−", "-").replace("–", "-").replace("π", "pi")
    t = re.sub(r"√\s*([\d.]+)", r"sqrt(\1)", t).replace("√", "sqrt")
    if not t:
        raise ValueError(f"{_contrae('a', nombre)} le falta un número: revisa que no sobren comas.")
    if not _CARACTERES.match(t):
        raise ValueError(f"en {nombre} hay un símbolo que no reconozco ({t!r}): usa solo números con punto "
                         "decimal, + - * / ^, paréntesis, sqrt(...) y pi.")
    if re.search(r"[\d.]\s+[\d.]", t):
        raise ValueError(f"en {nombre}, {t!r} tiene dos números separados solo por un espacio: "
                         "si son dos componentes, sepáralas con una coma.")
    if re.search(r"\.\d*\.", t):
        raise ValueError(f"en {nombre}, {t!r} tiene dos puntos decimales en un mismo número.")
    desconocidos = sorted(set(re.findall(r"[A-Za-z]\w*", _LITERAL.sub(" ", t))) - set(_PERMITIDOS))
    if desconocidos:
        raise ValueError(f"en {nombre} no reconozco {', '.join(desconocidos)}: escribe solo números, "
                         "fracciones como 2/3, sqrt(2) o pi.")
    try:
        expr = parse_expr(t, local_dict=dict(_PERMITIDOS), global_dict=dict(_GLOBAL),
                          transformations=_TRANSF, evaluate=False)
        # 9^9^9 no se evalúa (tardaría horas): se revisa cada potencia, de adentro hacia afuera
        grande = any(_potencia_excesiva(p) for p in sp.postorder_traversal(expr))
        valor = None
        if not grande:
            try:
                valor = sp.N(expr, 30, strict=True)
            except PrecisionExhausted:            # 4 - 4, pi - pi, 1/(1 - 1): se cancela; se evalúa exacto
                valor = sp.N(expr.doit(), 30)
    except Exception as err:
        raise ValueError(f"no pude leer {t!r} en {nombre}: usa * para multiplicar y revisa los paréntesis.") from err
    if grande:
        raise ValueError(f"en {nombre}, {t!r} tiene una potencia demasiado grande.")
    if not isinstance(valor, sp.Expr) or not valor.is_number or valor.has(sp.zoo, sp.nan, sp.oo, -sp.oo):
        raise ValueError(f"en {nombre}, {t!r} no da un número finito (¿dividiste entre 0?).")
    if sp.im(valor) != 0:
        raise ValueError(f"en {nombre}, {t!r} no es un número real (¿raíz de un negativo?).")
    valor = sp.re(valor)
    _revisa_rango(valor, nombre)
    return float(valor)


def numero(valor, nombre="el valor"):
    """Un número escrito como texto ('7', 'sqrt(49)', '2/3') o ya como número -> float."""
    if isinstance(valor, str):
        return _lee_numero(valor, nombre)
    if valor is None or isinstance(valor, (bool, np.bool_)):
        raise ValueError(f"{nombre} debe ser un número.")
    try:
        v = float(valor)
    except (TypeError, ValueError) as err:
        raise ValueError(f"{nombre} debe ser un número real.") from err
    if not math.isfinite(v):
        raise ValueError(f"{nombre} debe ser un número finito.")
    return v


def _sin_envoltura(t):
    """Quita un par de paréntesis o corchetes que encierra todo el texto: '(1, 2)' -> '1, 2'."""
    if len(t) < 2 or t[0] not in _CIERRE or t[-1] != _CIERRE[t[0]]:
        return t
    nivel = 0
    for i, ch in enumerate(t):
        nivel += (ch == t[0]) - (ch == _CIERRE[t[0]])
        if nivel == 0 and i < len(t) - 1:      # el primer paréntesis cierra antes del final: no envuelve todo
            return t
    return t[1:-1].strip()


def _partes(texto, nombre):
    t = texto.strip()
    while _sin_envoltura(t) != t:
        t = _sin_envoltura(t)
    if not t:
        raise ValueError(f"escribe {nombre} como números separados por comas, por ejemplo: 1, -2, 3")
    if t.count("(") != t.count(")"):
        raise ValueError(f"en {nombre} los paréntesis no cierran: revisa que cada ( tenga su ).")
    partes = [p.strip() for p in re.split(r"[,;]", t)] if re.search(r"[,;]", t) else t.split()
    if any(not p for p in partes):
        raise ValueError(f"en {nombre} sobra una coma o falta un número; escribe, por ejemplo: 1, -2, 3")
    return partes


def vector(texto, dim=None, nombre="el vector"):
    """'1, -2, 3', '(1,2,3)', '[1 2 3]', '2/3, sqrt(2)' o una lista de números -> arreglo de numpy (float).
    Con dim, exige esa cantidad de componentes."""
    if isinstance(texto, str):
        partes = _partes(texto, nombre)
    elif isinstance(texto, (list, tuple, np.ndarray)):
        partes = list(np.ravel(np.array(texto, dtype=object)))
    else:
        raise ValueError(f"escribe {nombre} como texto, por ejemplo: \"1, -2, 3\"")
    v = np.array([numero(p, nombre if isinstance(p, str) else f"cada componente {_contrae('de', nombre)}")
                  for p in partes], dtype=float)
    if dim is not None and len(v) != dim:
        raise ValueError(f"{nombre} debe tener {dim} componentes y tiene {len(v)}.")
    if len(v) < 2:
        raise ValueError(f"{nombre} necesita al menos dos componentes separadas por comas, por ejemplo: 1, -2, 3")
    return v


def cifras(valor, n=4):
    """Redondea (no trunca) a n cifras significativas con la regla escolar (5 sube) y conserva los ceros finales.
    Con exponente decimal mayor que 15 o menor que -4 usa notación científica: 1.234×10^-7."""
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
    e = q.adjusted()
    if e > 15 or e < -4:                      # muy grande o muy pequeño: mantisa × 10^e
        return f"{format(q.scaleb(-e), 'f')}×10^{e}"
    if e >= n:                                # enteros grandes: sin notación científica
        return f"{int(q):d}"
    return format(q, "f")


def fmt(v, n=4):
    """Vector como texto '(a, b, c)', cada componente redondeada a n cifras significativas."""
    return "(" + ", ".join(cifras(c, n) for c in np.ravel(v)) + ")"


def cerca(p, q, tol=1e-9):
    return abs(float(p) - float(q)) < tol


def _rechaza(fn):
    """True si fn() lanza ValueError (para probar que una calculadora rechaza una entrada)."""
    try:
        fn()
    except ValueError:
        return True
    return False


def en_unidad(rad, unidad):
    """Ángulo en radianes -> en la unidad pedida ('grados' o 'radianes')."""
    if unidad not in ("grados", "radianes"):
        raise ValueError("La unidad debe ser 'grados' o 'radianes'.")
    return math.degrees(rad) if unidad == "grados" else float(rad)


def de_unidad(valor, unidad):
    """Ángulo en 'grados' o 'radianes' -> en radianes."""
    if unidad not in ("grados", "radianes"):
        raise ValueError("La unidad debe ser 'grados' o 'radianes'.")
    return math.radians(valor) if unidad == "grados" else float(valor)


def _sufijo(unidad):
    return "°" if unidad == "grados" else " rad"


# --- cálculo interno (sin validar: lo usan las funciones públicas después de validar) ---
def _norma(v):
    return math.hypot(*v)                     # hypot no se desborda al elevar al cuadrado


def _punto(a, b):
    """a · b con la cancelación por redondeo limpiada: 1e-17 que debería ser 0 se vuelve 0."""
    productos = [float(x) * float(y) for x, y in zip(a, b)]
    s = math.fsum(productos)
    return 0.0 if abs(s) <= 1e-14 * math.fsum(abs(p) for p in productos) else s


def _limpia(valores, escalas):
    """Pone en 0 cada componente que solo es ruido de redondeo frente a su escala."""
    v = np.array(valores, dtype=float)
    v[np.abs(v) <= 1e-14 * np.asarray(escalas, dtype=float)] = 0.0
    return v + 0.0                            # + 0.0 cambia -0.0 por 0.0


def _escalado(v):
    """v multiplicado por una potencia de 2 (exacto) para que su mayor componente quede entre 0.5 y 1."""
    m = float(np.max(np.abs(v)))
    return v if m == 0 else np.ldexp(v, -math.frexp(m)[1])


def _par(a, b):
    """Lee a y b y exige que tengan la misma cantidad de componentes."""
    a, b = vector(a, nombre="a"), vector(b, nombre="b")
    if len(a) != len(b):
        raise ValueError(f"a y b deben tener la misma cantidad de componentes; tienen {len(a)} y {len(b)}.")
    return a, b


def _ejes_iguales(ax, puntos, margen=0.15):
    """Misma escala en todos los ejes (2D o 3D) para que los ángulos no se deformen."""
    P = np.atleast_2d(np.array(puntos, dtype=float))
    bajo, alto = P.min(axis=0), P.max(axis=0)
    centro = (bajo + alto) / 2
    radio = (alto - bajo).max() / 2 * (1 + margen)
    if not radio > 0:
        radio = max(1.0, 0.1 * np.abs(centro).max())
    lims = [(c - radio, c + radio) for c in centro]
    ax.set_xlim(*lims[0]); ax.set_ylim(*lims[1])
    if P.shape[1] == 3:
        ax.set_zlim(*lims[2]); ax.set_box_aspect((1, 1, 1))
    else:
        ax.set_aspect("equal")


def _vista(ax, vectores=(), normal=None):
    """Elige la vista 3D (elevación de 15° a 40°) en la que cada vector queda a 30° o más de la línea de
    visión y, si hay normal, su plano se ve inclinado (de 30° a 60° respecto a verlo de frente). Entre las
    vistas que lo cumplen, la más cercana a la de matplotlib por omisión (azimut -60°, elevación 30°)."""
    def unitario_o_nada(v):
        v = np.asarray(v, dtype=float)
        return v / _norma(v) if np.any(v) and np.all(np.isfinite(v)) else None
    dirs = [d for d in map(unitario_o_nada, vectores) if d is not None]
    nn = None if normal is None else unitario_o_nada(normal)
    mejor = None
    for elev in range(15, 41, 5):
        for azim in range(-180, 180, 5):
            e, a = math.radians(elev), math.radians(azim)
            ojo = np.array([math.cos(e) * math.cos(a), math.cos(e) * math.sin(a), math.sin(e)])
            senos = [math.sqrt(max(0.0, 1 - (ojo @ d) ** 2)) for d in dirs]
            if nn is not None:
                senos += [math.sqrt(max(0.0, 1 - (ojo @ nn) ** 2)), abs(ojo @ nn)]
            lejania = abs(elev - 30) + abs((azim + 60 + 180) % 360 - 180) / 4
            clave = (min(round(min(senos, default=1.0), 2), 0.5), -lejania)   # sen 30° = 0.5 basta
            if mejor is None or clave > mejor[0]:
                mejor = (clave, elev, azim)
    ax.view_init(elev=mejor[1], azim=mejor[2])
