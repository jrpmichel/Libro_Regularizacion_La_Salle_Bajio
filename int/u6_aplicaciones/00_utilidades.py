# ID: INT-U6-NB00
# Notebook: int/u6_aplicaciones.ipynb · utilidades compartidas
# Repositorio: int/u6_aplicaciones/00_utilidades.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# # Int. U6 · Aplicaciones a ingeniería
# **Libro de regularización en matemáticas para ingeniería · ULSA Bajío**
#
# Este notebook acompaña la unidad 6 del bloque de cálculo integral. Tiene una sección por subtema, en el mismo orden que el libro:
#
# - 6.1 Área entre curvas
# - 6.2 Volúmenes de sólidos de revolución
# - 6.3 Longitud de arco
# - 6.4 Centroides y momentos de inercia de área
# - 6.5 Momentos de masa
# - 6.6 De la aceleración a la posición
# - 6.7 Trabajo, energía y eficiencia
# - 6.8 Fuerza hidrostática y caudal
# - 6.9 Carga, energía y valor eficaz de una señal
# - 6.10 Costo acumulado
#
# **Cómo usarlo.** Ejecuta las celdas en orden (Entorno de ejecución → Ejecutar todas). Cada sección tiene una celda de teoría con un ejemplo exacto, una calculadora, sus casos de prueba y una celda **Contrasta**. Antes de ejecutar la celda Contrasta, resuelve a mano el caso que se indica y escribe tu resultado. La referencia de la calculadora aparece solo después de que escribes el tuyo.
#
# **Etiqueta del repositorio.** `INT-U6` (bloque, unidad). Las celdas de código de cada sección vienen de `int/u6_aplicaciones/NN_*.py`; los fragmentos que el libro imprime, de `int/u6_aplicaciones/fragmentos/`.
#
# **Convenciones.** Punto decimal. Los valores se *redondean* (no se truncan) a las cifras significativas que declara cada calculadora. Los ángulos de `sin`, `cos` y `tan` van en **radianes**. Cada resultado lleva su unidad: las secciones de física usan el SI (con prefijos como mA o µA cuando el dato viene así) y la sección 6.10 usa pesos. Para escribir funciones usa `x` como variable o la letra natural de la sección (`t`, `y`, `h`, `r` o `q`), `*` para multiplicar, `**` o `^` para potencias, `sqrt(...)`, `exp(...)`, `log(...)` o `ln(...)` (natural), `abs(...)`, `sen` o `sin`, `cosh`, `pi` o `π`, `√x`, números como `2e5` y `escalon(t - 2)` para un escalón unitario que vale 0 antes de $t=2$ y 1 después.

# %% [markdown]
# ## Versiones e imports

# %%
# Versiones mínimas probadas. En Colab ya vienen instaladas; esta celda solo lo asegura.
%pip install -q "numpy>=1.26" "sympy>=1.12" "matplotlib>=3.8" "ipywidgets>=8" "scipy>=1.11"

import sys, math, re, signal, warnings
from decimal import Decimal, ROUND_HALF_UP
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sympy as sp
import ipywidgets as widgets
from scipy.integrate import quad, IntegrationWarning
from scipy.optimize import brentq
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                        implicit_multiplication_application, convert_xor, rationalize)

import scipy
print("Python", sys.version.split()[0], "| numpy", np.__version__, "| scipy", scipy.__version__, "| sympy", sp.__version__,
      "| matplotlib", matplotlib.__version__, "| ipywidgets", widgets.__version__)

# %%
# Símbolos y utilidades compartidas por todo el notebook
x = sp.symbols("x", real=True)
CIFRAS = 6                                          # cifras significativas de las calculadoras

_FUNCIONES = {"sqrt": sp.sqrt, "abs": sp.Abs, "exp": sp.exp, "log": sp.log, "ln": sp.log,
              "pi": sp.pi, "sin": sp.sin, "sen": sp.sin, "cos": sp.cos, "tan": sp.tan,
              "atan": sp.atan, "arctan": sp.atan, "asin": sp.asin, "arcsen": sp.asin, "arcsin": sp.asin,
              "sinh": sp.sinh, "senh": sp.sinh, "cosh": sp.cosh, "escalon": sp.Heaviside, "Heaviside": sp.Heaviside,
              "E": sp.E, "e": sp.E}
_GLOBAL = {"Integer": sp.Integer, "Float": sp.Float, "Rational": sp.Rational,
           "Symbol": sp.Symbol, "Function": sp.Function}
# rationalize: 0.1 se lee como 1/10, para que los resultados salgan exactos
_TRANSF = standard_transformations + (implicit_multiplication_application, convert_xor, rationalize)
_OTRAS_VARIABLES = ("t", "y", "h", "r", "q", "s")     # se leen como x si son la única variable
_LITERAL = r"(?<![A-Za-z_\d.])(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?"   # 2, 0.5, .5, 2e5, 1.5e-3


def parsear(texto):
    """Convierte un texto en una expresión de sympy en x o explica qué salió mal.
    Acepta sen y ln, la letra e como número de Euler, π y √ (2πx, √x, √(x + 1)), |...| como valor absoluto,
    sin^2(x), números como 2e5, escalon(...) como escalón unitario y otra letra (t, y, h, r, q o s) en lugar de x
    si es la única variable."""
    if not isinstance(texto, str) or not texto.strip():
        raise ValueError("Escribe una función, por ejemplo: x**2")
    if "," in texto:
        raise ValueError("Usa punto decimal (1.5, no 1,5); las funciones de este notebook no llevan comas.")
    texto = texto.replace("·", "*").replace("−", "-").replace("escalón", "escalon")
    texto = texto.replace("π", " pi ")                                     # 2πx -> 2 pi x
    texto = re.sub(r"√\s*(\d+\.?\d*|[A-Za-z_]\w*|\([^()]*\))", r" sqrt(\1) ", texto)   # √x, √2, √(x + 1)
    texto = re.sub(r"\|([^|]+)\|", r"abs(\1)", texto)                      # |x - 1| -> abs(x - 1)
    texto = re.sub(r"\b(sin|sen|cos|tan)\s*\^\s*(\d+)\s*\(([^()]*)\)", r"(\1(\3))**\2", texto)   # sin^2(x)
    if "√" in texto:
        raise ValueError("Escribe la raíz como sqrt(...), por ejemplo sqrt(x + 1).")
    def _nombres(t):
        return set(re.findall(r"[A-Za-z_]\w*", re.sub(_LITERAL, " ", t)))   # sin los números: 2e5 no es un nombre
    nombres = _nombres(texto)
    otras = [v for v in _OTRAS_VARIABLES if v in nombres]
    if "x" not in nombres and len(otras) == 1:
        texto = re.sub(rf"(?<![A-Za-z_]){otras[0]}(?![A-Za-z_0-9])", "x", texto)
        nombres = _nombres(texto)
    local = {**_FUNCIONES, "x": x}
    desconocidos = sorted(nombres - set(local))
    if desconocidos:
        raise ValueError(f"No reconozco {desconocidos}. Usa una sola letra como variable (x, t, y, h, r o q), números "
                         "en lugar de letras para las constantes y, si hace falta, sqrt, exp, log o ln, abs, sin o sen, "
                         "cos, tan, cosh, escalon o pi.")
    try:
        expr = parse_expr(texto, local_dict=local, global_dict=dict(_GLOBAL), transformations=_TRANSF)
    except Exception as err:
        raise ValueError("No pude leer la función. Usa * para multiplicar, ** o ^ para potencias, punto decimal "
                         "y sqrt(...) para raíces.") from err
    if not isinstance(expr, sp.Expr) or expr.atoms(sp.core.function.AppliedUndef):
        raise ValueError("Escribe una expresión sin el signo = y solo con las funciones permitidas.")
    if expr.has(sp.oo, -sp.oo, sp.zoo, sp.nan):
        raise ValueError("La función no puede contener infinito.")
    if expr.free_symbols - {x}:
        raise ValueError("Usa una sola variable.")
    return expr


def lista_numeros(texto, nombre="la lista"):
    """'12, 2, 11' -> [12.0, 2.0, 11.0], o un mensaje claro si algo no es número."""
    partes = [p for p in re.split(r"[,\s;]+", str(texto).strip()) if p]
    if not partes:
        raise ValueError(f"Escribe {nombre} como números separados por comas, por ejemplo: 12, 2, 11")
    try:
        return [float(p) for p in partes]
    except ValueError as err:
        raise ValueError(f"En {nombre} hay algo que no es un número: usa punto decimal y comas entre valores.") from err


def cifras(valor, n=4):
    """Redondea (no trunca) a n cifras significativas con la regla escolar (5 sube) y conserva los ceros finales.
    Fuera de 1e-4 a 1e7 usa notación científica."""
    v = float(valor)
    if math.isnan(v):
        return "no definido"
    if math.isinf(v):
        return "+∞" if v > 0 else "-∞"
    if v == 0:
        return "0"
    d = Decimal(repr(v))
    if abs(v) >= 1e7 or abs(v) < 1e-4:
        m = (d.scaleb(-d.adjusted())).quantize(Decimal(1).scaleb(-n + 1), rounding=ROUND_HALF_UP)
        e = d.adjusted()
        if abs(m) >= 10:                     # 9.996e5 con 3 cifras pasa a 1.00e6
            m, e = (m / 10).quantize(Decimal(1).scaleb(-n + 1), rounding=ROUND_HALF_UP), e + 1
        return f"{m}e{e}"
    q = d.quantize(Decimal(1).scaleb(d.adjusted() - n + 1), rounding=ROUND_HALF_UP)
    if q.adjusted() != d.adjusted():          # 9.996 con 3 cifras pasa a 10.0: se reajusta la posición
        q = d.quantize(Decimal(1).scaleb(q.adjusted() - n + 1), rounding=ROUND_HALF_UP)
    if q.adjusted() >= n:                     # enteros grandes: sin notación científica
        return f"{int(q):d}"
    return format(q, "f")


def cerca(p, q, tol=1e-9):
    return abs(float(p) - float(q)) < tol


def coincide(mi, ref, n=3):
    """True si mi y ref coinciden a n cifras significativas (para las celdas Contrasta)."""
    mi, ref = float(mi), float(ref)
    return abs(mi - ref) <= 5 * 10**(-n) * max(abs(ref), 1e-12) or abs(mi - ref) < 1e-12


def a_numero(valor):
    """Número escrito por el alumno: 1.57, 'pi/2' o '3/2'. ValueError con un mensaje claro si no es un real finito."""
    if isinstance(valor, bool):
        raise ValueError("Escribe un número, no True ni False.")
    if isinstance(valor, (int, float)):
        if not math.isfinite(valor):
            raise ValueError("Escribe un número finito.")
        return float(valor)
    texto = str(valor).strip()
    if not texto:
        raise ValueError("Escribe tu resultado entre comillas, por ejemplo \"1.5\" o \"4/3\".")
    if re.fullmatch(r"[+-]?\d*,\d+", texto):
        raise ValueError(f"Usa punto decimal: escribe {texto.replace(',', '.')} en lugar de {texto}.")
    try:
        v = parsear(texto)
    except ValueError as err:
        if "infinito" in str(err):
            raise ValueError("Eso no es un número finito (¿dividiste entre cero?).") from None
        raise
    if v.free_symbols:
        raise ValueError("Escribe un número, sin x.")
    try:
        z = complex(sp.N(v))
    except (TypeError, ValueError):
        raise ValueError("Eso no es un número real.") from None
    if abs(z.imag) > 1e-12 or not math.isfinite(z.real):
        raise ValueError("Eso no es un número real finito.")
    return z.real


def _rechaza(fn):
    """True si fn() lanza ValueError (para probar que una calculadora rechaza una entrada)."""
    try:
        fn()
    except ValueError:
        return True
    return False


class _Tarda(BaseException):                 # BaseException: los "except Exception" internos de sympy no la atrapan
    pass


def con_limite(fn, segundos=15):
    """Ejecuta fn() con un límite de tiempo (en Colab); si se pasa, lanza ValueError con un mensaje claro."""
    if not hasattr(signal, "SIGALRM"):        # Windows: sin límite
        return fn()
    def _alarma(*_):
        raise _Tarda()
    anterior = signal.signal(signal.SIGALRM, _alarma)
    signal.alarm(segundos)
    try:
        return fn()
    except _Tarda:
        raise ValueError(f"sympy tardó más de {segundos} s con esta función; prueba una más sencilla.") from None
    finally:
        signal.alarm(0)
        signal.signal(signal.SIGALRM, anterior)


# Errores que una calculadora convierte en un mensaje para el alumno
ERRORES = (ValueError, TypeError, ArithmeticError)


def explica(err):
    """Mensaje en español para un error de cálculo."""
    if isinstance(err, ValueError):
        return str(err)
    if isinstance(err, OverflowError):
        return "los números crecen demasiado (desbordamiento): usa un intervalo más corto."
    if isinstance(err, ZeroDivisionError):
        return "apareció una división entre cero: revisa los datos."
    return "no pude convertir el resultado en un número real; prueba con otra función o con otro intervalo."


def numerica(f):
    """Versión numpy de f que siempre devuelve un arreglo de floats del tamaño de la entrada (nan donde no es real)."""
    g = sp.lambdify(x, sp.sympify(f), "numpy")
    def fn(t):
        t = np.asarray(t, dtype=float)
        with np.errstate(all="ignore"):
            y = np.asarray(g(t)) * np.ones_like(t)
        if np.iscomplexobj(y):
            y = np.where(np.abs(y.imag) < 1e-12, y.real, np.nan)
        return np.asarray(y, dtype=float)
    return fn


def para_graficar(f, xs):
    """Valores de f en xs, con nan donde no son finitos (matplotlib deja un hueco en lugar de fallar)."""
    ys = numerica(f)(xs)
    return np.where(np.isfinite(ys), ys, np.nan)


def exacto(v, nombre="el valor"):
    """Número de sympy, exacto si se puede (0.1 -> 1/10, 'pi/2' -> pi/2); ValueError si no es un real finito."""
    if isinstance(v, bool):
        raise ValueError(f"{nombre} debe ser un número.")
    if isinstance(v, sp.Basic):
        e = v
    elif isinstance(v, (int, np.integer)):
        e = sp.Integer(int(v))
    elif isinstance(v, (float, np.floating)):
        if not math.isfinite(float(v)):
            raise ValueError(f"{nombre} debe ser un número finito.")
        e = sp.Rational(repr(float(v)))
    else:
        try:
            e = parsear(str(v))
        except ValueError:
            raise ValueError(f"Escribe {nombre} como número: 3, 0.5, 1/2 o pi/2.") from None
    if getattr(e, "free_symbols", None):
        raise ValueError(f"{nombre} es un número: no lleva variable.")
    try:
        z = complex(sp.N(e))
    except (TypeError, ValueError):
        raise ValueError(f"{nombre} debe ser un número real.") from None
    if abs(z.imag) > 0 or not math.isfinite(z.real):
        raise ValueError(f"{nombre} debe ser un número real finito.")
    return e


def simplifica(v, segundos=5):
    """sympy.simplify con límite de tiempo; si tarda, deja el valor como está."""
    if not isinstance(v, sp.Basic) or isinstance(v, sp.Float):
        return v
    try:
        s_ = con_limite(lambda: sp.simplify(v), segundos)
    except ValueError:
        return v
    e_ = sp.expand(s_)                           # a veces la forma expandida es más corta: 1/100 - exp(-5)/100
    return e_ if len(str(e_)) < len(str(s_)) else s_


def _es_decimal_finito(v):
    """True si el racional v tiene expansión decimal finita (denominador con solo 2 y 5)."""
    q = int(v.q)
    for p_ in (2, 5):
        while q % p_ == 0:
            q //= p_
    return q == 1


def muestra(v, n=CIFRAS):
    """El valor redondeado a n cifras significativas. Si v es exacto y no es un decimal finito (16/3, 2*sqrt(2), pi),
    antepone su forma exacta: '16/3 ≈ 5.33333'. Un entero corto se escribe tal cual."""
    if isinstance(v, sp.Basic) and v.atoms(sp.Float):     # ya es decimal (viene de quad o de un cruce numérico)
        return cifras(v, n)
    if isinstance(v, sp.Integer) and len(str(abs(int(v)))) <= n:
        return str(v)
    if isinstance(v, sp.Rational) and _es_decimal_finito(v):
        return cifras(v, n)
    if isinstance(v, sp.Basic) and not isinstance(v, (sp.Float, sp.Integer)) and len(str(v)) < 40:
        return f"{v} ≈ {cifras(v, n)}"
    return cifras(v, n)


def malla_interior(a, b, n=401):
    """n puntos dentro de (a, b), sin los extremos (ahí puede haber una singularidad integrable)."""
    return np.linspace(float(a), float(b), n + 2)[1:-1]


def revisa_intervalo(f, a, b, nombre="f", var="x"):
    """ValueError si f no es real y continua dentro de (a, b). Decide sympy; con escalones, o si sympy no puede,
    se revisa en una malla de puntos. Los extremos se dejan libres: ahí la integral puede ser impropia y convergente.
    var es la letra que ve el alumno en los mensajes (t, h, r, q)."""
    lo, hi = sorted((a, b), key=float)
    f = sp.sympify(f)
    if not f.has(x):
        return
    if not f.has(sp.Heaviside, sp.Piecewise):
        try:
            dom = con_limite(lambda: sp.calculus.util.continuous_domain(f, x, sp.Interval.open(lo, hi)), 5)
        except (ValueError, NotImplementedError, TypeError):
            dom = None
        if dom is not None and dom != sp.Interval.open(lo, hi):
            fuera = sp.Interval.open(lo, hi) - dom
            if isinstance(fuera, sp.FiniteSet) and len(fuera) <= 4:
                donde = ", ".join(punto(p) for p in sorted(fuera, key=float))
                raise ValueError(f"{nombre} no es continua en {var} = {donde}: ahí hay una asíntota o un salto. "
                                 "Separa el intervalo o elige otro.")
    ys = numerica(f)(malla_interior(lo, hi))
    if np.any(np.isnan(ys)):
        raiz_impar = any(p.exp.is_Rational and p.exp.q % 2 == 1 and p.exp.q > 1 for p in f.atoms(sp.Pow))
        raise ValueError(f"{nombre} no está definida (o no es real) en todo el intervalo: revisa el dominio."
                         + (" Para sympy, una potencia como x^(2/3) o x^(1/3) no es real si la base es negativa; "
                            "usa un intervalo donde la base sea positiva." if raiz_impar else ""))
    if np.any(np.isinf(ys)):
        raise ValueError(f"{nombre} crece demasiado dentro del intervalo: sus valores pasan de 1e308 (desbordamiento). "
                         "Usa un intervalo más corto.")


# Integrales definidas: exactas con sympy si se puede; si no, numéricas con quad
_RARAS = (sp.Integral, sp.exp_polar, sp.hyper, sp.meijerg, sp.lowergamma, sp.uppergamma, sp.Piecewise,
          sp.polar_lift, sp.Min, sp.Max)


def _valor_exacto(valor):
    """float real del resultado de sympy, o None si no sirve: no es un número, es complejo, usa funciones que sympy
    no compara bien (exp_polar, hyper, ...) o tiene condiciones sin decidir."""
    if not isinstance(valor, sp.Expr) or valor.free_symbols or valor.has(*_RARAS) or valor.has(sp.I):
        return None
    if valor.atoms(sp.core.relational.Relational):
        return None
    try:
        z = complex(con_limite(lambda: sp.N(valor), 5))
    except (ValueError, TypeError, ArithmeticError):
        return None
    if not (math.isfinite(z.real) and math.isfinite(z.imag)) or abs(z.imag) > 1e-12 * max(1.0, abs(z.real)):
        return None
    return z.real


def _quad(fn, a, b):
    """∫_a^b por cuadratura, en tramos de longitud 50 como máximo (las oscilaciones largas no confunden a quad)."""
    tramos = min(2000, max(1, int(math.ceil(abs(b - a) / 50))))
    cortes = np.linspace(a, b, tramos + 1)
    total = 0.0
    with warnings.catch_warnings():
        warnings.simplefilter("error", IntegrationWarning)
        try:
            for p, q in zip(cortes, cortes[1:]):
                total += quad(lambda s: float(fn(s)), p, q, limit=200)[0]
        except (IntegrationWarning, ZeroDivisionError, OverflowError):
            raise ValueError("quad no alcanzó la precisión pedida: la integral puede diverger (si la función crece sin "
                             "límite cerca de algún punto) o la función oscila demasiado. Prueba con un intervalo más "
                             "corto.") from None
    if not math.isfinite(total):
        raise ValueError("la integral numérica no dio un número finito: revisa el dominio de la función.")
    return total


def integra(f, a, b, segundos=8, solo_numerica=False):
    """∫_a^b f(x) dx. Devuelve (valor, "exacta") si sympy la resuelve en el tiempo dado con una expresión útil;
    si no, (Float, "numérica") con quad. ValueError si diverge o si f no es real en el intervalo."""
    f = sp.sympify(f)
    a, b = exacto(a, "el extremo inferior"), exacto(b, "el extremo superior")
    if f.has(sp.Heaviside):                     # sympy integra mejor los escalones como funciones por tramos
        f = f.rewrite(sp.Piecewise)
    valor = None
    if not solo_numerica:
        try:
            valor = con_limite(lambda: sp.integrate(f, (x, a, b)), segundos)
        except (ValueError, NotImplementedError, TypeError, ArithmeticError):
            valor = None
    if isinstance(valor, sp.Expr) and (valor.has(sp.oo, -sp.oo, sp.zoo) and not valor.has(sp.nan)):
        raise ValueError("la integral diverge: la función crece sin límite en el intervalo.")
    z = _valor_exacto(valor) if valor is not None else None
    if z is not None:
        if valor.atoms(sp.Float) or any(abs(r.q) > 10**6 for r in valor.atoms(sp.Rational)):
            return sp.Float(z, 15), "exacta"     # extremos decimales largos (un cruce numérico): el valor en decimal
        return simplifica(valor), "exacta"
    return sp.Float(_quad(numerica(f), float(a), float(b)), 15), "numérica"


class Integrador:
    """Integra varias funciones en una misma llamada. En cuanto sympy falla o tarda una vez, las demás integrales
    van directo a quad; así una calculadora no espera el límite de tiempo en cada integral."""
    def __init__(self, segundos=8):
        self.segundos, self.solo_numerica, self.metodos = segundos, False, set()

    def __call__(self, f, a, b):
        v, m = integra(f, a, b, self.segundos, self.solo_numerica)
        self.metodos.add(m)
        self.solo_numerica = self.solo_numerica or m == "numérica"
        return v

    @property
    def metodo(self):
        return "numérica" if "numérica" in self.metodos else "exacta"


def cruces(h, a, b):
    """Puntos de (a, b) donde h cambia de signo: exactos con solveset; si no, en una malla y con brentq.
    Devuelve (lista, "exactos" o "numéricos")."""
    try:
        sol = con_limite(lambda: sp.solveset(h, x, sp.Interval.open(a, b)), 8)
    except (ValueError, NotImplementedError, TypeError):
        sol = None
    if sol == sp.EmptySet:
        return [], "exactos"
    if isinstance(sol, sp.FiniteSet) and len(sol) <= 50 and all(c.is_real for c in sol):
        return sorted(sol, key=float), "exactos"
    fn = numerica(h)
    xs = malla_interior(a, b, 4001)
    signo = np.sign(fn(xs))
    raices, previo = [], None                       # previo: índice del último punto con signo distinto de 0
    for k, s_ in enumerate(signo):
        if s_ == 0 or np.isnan(s_):
            continue
        if previo is not None and s_ != signo[previo]:
            if k == previo + 1:
                raices.append(brentq(lambda t: float(fn(t)), xs[previo], xs[k], xtol=1e-14))
            else:                                   # h = 0 justo en la malla
                raices.append(xs[(previo + k) // 2])
        previo = k
    return [sp.Float(r, 15) for r in raices], "numéricos"


COMO = {"exacta": "cálculo exacto con sympy", "numérica": "cálculo numérico con quad"}


def punto(v, n=CIFRAS):
    """Un extremo o un cruce legible: exacto si es corto (pi/4, 3/2); si no, redondeado a n cifras."""
    if isinstance(v, sp.Basic) and not isinstance(v, sp.Float) and len(str(v)) <= 12:
        return str(v)
    return cifras(v, n)


def rotula(ax, xv, yv, texto, **kw):
    """Etiqueta directa junto a un punto de una curva (legible en blanco y negro, sin leyenda)."""
    kw = {"xytext": (4, 0), "textcoords": "offset points", "va": "center", "fontsize": 9, **kw}
    ax.annotate(texto, (float(xv), float(yv)), **kw)
