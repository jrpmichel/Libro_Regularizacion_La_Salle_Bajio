# ID: INT-U5-NB00
# Notebook: int/u5_tecnicas.ipynb · utilidades compartidas
# Repositorio: int/u5_tecnicas/00_utilidades.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# # Int. U5 · Técnicas de integración
# **Libro de regularización en matemáticas para ingeniería · ULSA Bajío**
#
# Este notebook acompaña la unidad 5 del bloque de cálculo integral. Tiene una sección por subtema, en el mismo orden que el libro:
#
# 1. Sustitución
# 2. Integración por partes
# 3. Fracciones parciales
# 4. Potencias de funciones trigonométricas y sustitución trigonométrica
# 5. Integrales impropias
#
# **Cómo usarlo.** Ejecuta las celdas en orden (Entorno de ejecución → Ejecutar todas). Cada sección termina con una celda **Contrasta**: antes de ejecutarla, resuelve a mano el caso que se indica y escribe tu resultado. La referencia de la calculadora aparece solo después de que escribes el tuyo.
#
# **Etiqueta del repositorio.** `INT-U5` (bloque, unidad). Las celdas de código de cada sección vienen de `int/u5_tecnicas/NN_*.py`; los fragmentos que el libro imprime, de `int/u5_tecnicas/fragmentos/`.
#
# **Convenciones.** Punto decimal. Los valores se *redondean* (no se truncan) a las cifras significativas que elijas. Para escribir funciones usa `x` como variable, `*` para multiplicar, `**` o `^` para potencias, `sqrt(...)`, `exp(...)`, `log(...)` (natural), `abs(...)`, `atan(...)` y `asin(...)`; también `pi`, `sen` y `ln`. Los ángulos de `sin`, `cos` y `tan` van en **radianes**. Toda antiderivada se muestra con su constante `+ C`. Las calculadoras de antiderivadas derivan el resultado para comprobarlo y lo comparan con la antiderivada de `sympy`, porque dos respuestas correctas pueden diferir en una constante.

# %% [markdown]
# ## 0. Versiones e imports

# %%
# Versiones mínimas probadas. En Colab ya vienen instaladas; esta celda solo lo asegura.
%pip install -q "numpy>=1.26" "sympy>=1.12" "matplotlib>=3.8" "ipywidgets>=8" "scipy>=1.11"

import sys, math, re, signal
from decimal import Decimal, ROUND_HALF_UP
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import sympy as sp
import ipywidgets as widgets
from IPython.display import display
from scipy.integrate import quad
from scipy.optimize import brentq
from sympy.integrals import manualintegrate as _mi
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                        implicit_multiplication_application, convert_xor, rationalize)

import scipy
print("Python", sys.version.split()[0], "| numpy", np.__version__, "| scipy", scipy.__version__, "| sympy", sp.__version__,
      "| matplotlib", matplotlib.__version__, "| ipywidgets", widgets.__version__)

# %%
# Símbolos y utilidades compartidas por todo el notebook
x = sp.symbols("x", real=True)
U = sp.symbols("u", real=True)                     # variable de las sustituciones

_FUNCIONES = {"sqrt": sp.sqrt, "abs": sp.Abs, "exp": sp.exp, "log": sp.log, "ln": sp.log,
              "pi": sp.pi, "sin": sp.sin, "sen": sp.sin, "cos": sp.cos, "tan": sp.tan, "sec": sp.sec,
              "atan": sp.atan, "arctan": sp.atan, "asin": sp.asin, "arcsen": sp.asin, "arcsin": sp.asin,
              "acos": sp.acos, "arccos": sp.acos, "sinh": sp.sinh, "cosh": sp.cosh, "E": sp.E, "e": sp.E}
_GLOBAL = {"Integer": sp.Integer, "Float": sp.Float, "Rational": sp.Rational,
           "Symbol": sp.Symbol, "Function": sp.Function}
# rationalize: 0.1 se lee como 1/10, para que los resultados salgan exactos
_TRANSF = standard_transformations + (implicit_multiplication_application, convert_xor, rationalize)
_OTRAS_VARIABLES = ("t", "w", "q", "s")             # se leen como x si son la única variable (u se reserva)
_C = sp.Symbol("C")
_c = sp.Symbol("c")                               # también se acepta + c minúscula


def parsear(texto, variables=("x",)):
    """Convierte un texto en una expresión de sympy o explica qué salió mal.
    Acepta sen y ln, la letra e como número de Euler, π y √, |...| como valor absoluto, sin^2(x),
    una constante + C o + c (se descarta) y otra letra (t, w, q o s) en lugar de x si es la única variable."""
    if not isinstance(texto, str) or not texto.strip():
        raise ValueError("Escribe una función, por ejemplo: x**2")
    texto = texto.replace("π", "pi").replace("√", "sqrt").replace("·", "*").replace("−", "-")
    texto = re.sub(r"\|([^|]+)\|", r"abs(\1)", texto)                      # |x - 1| -> abs(x - 1)
    texto = re.sub(r"\b(sin|sen|cos|tan)\s*\^\s*(\d+)\s*\(([^()]*)\)", r"(\1(\3))**\2", texto)   # sin^2(x)
    nombres = set(re.findall(r"[A-Za-z_]\w*", texto))
    otras = [v for v in _OTRAS_VARIABLES if v in nombres]
    if "x" not in nombres and len(otras) == 1:
        texto = re.sub(rf"(?<![A-Za-z_]){otras[0]}(?![A-Za-z_0-9])", "x", texto)
        nombres = set(re.findall(r"[A-Za-z_]\w*", texto))
    local = {**_FUNCIONES, "x": x, "C": _C, "c": _c}
    desconocidos = sorted(nombres - set(local))
    if desconocidos:
        raise ValueError(f"No reconozco {desconocidos}. Usa x (o t) como variable, números en lugar de letras para las "
                         "constantes y, si hace falta, sqrt, exp, log o ln, abs, sin o sen, cos, tan, atan, asin o pi.")
    try:
        expr = parse_expr(texto, local_dict=local, global_dict=dict(_GLOBAL), transformations=_TRANSF)
    except Exception as err:
        raise ValueError("No pude leer la función. Usa * para multiplicar, ** o ^ para potencias, punto decimal "
                         "y sqrt(...) para raíces.") from err
    if not isinstance(expr, sp.Expr) or expr.atoms(sp.core.function.AppliedUndef):
        raise ValueError("Escribe una expresión sin el signo = y solo con las funciones permitidas.")
    expr = expr.subs({_C: 0, _c: 0})              # la constante de integración no cambia la derivada
    if expr.has(sp.oo, -sp.oo, sp.zoo, sp.nan):
        raise ValueError("La función no puede contener infinito; los extremos infinitos se escriben aparte.")
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
    """Número escrito por el alumno: 1.57, 'pi/2' o '3/2'."""
    if isinstance(valor, (int, float)):
        return float(valor)
    v = parsear(str(valor))
    if v.free_symbols:
        raise ValueError("Escribe un número, sin x.")
    return float(sp.N(v))


def _rechaza(fn):
    """True si fn() lanza ValueError (para probar que una calculadora rechaza una entrada)."""
    try:
        fn()
    except ValueError:
        return True
    return False


class _Tarda(Exception):
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


# Antiderivadas (las mismas herramientas de las Unidades 3 y 4, ampliadas para esta unidad)
_PERMITIDAS = (sp.exp, sp.log, sp.sin, sp.cos, sp.tan, sp.sec, sp.Abs, sp.atan, sp.asin, sp.acos)
_PUNTOS = [-8.9, -5.7, -4.3, -2.9, -2.1, -1.3, -0.7, -0.37, 0.13, 0.37, 0.53, 0.7, 0.91, 1.3, 2.1, 2.9, 4.3, 5.7, 8.9]


def _siempre_positivo(a):
    """True si a > 0 para todo x real (sympy lo sabe, o es un polinomio sin raíces reales y coeficiente principal > 0)."""
    if a.is_positive:
        return True
    if a.has(x) and a.is_polynomial(x):
        P = sp.Poly(a, x)
        return P.count_roots() == 0 and P.LC() > 0
    return False


def con_valor_absoluto(F, conservar=()):
    """ln(u) -> ln|u| salvo que u sea siempre positivo: la antiderivada de u'/u vale donde u > 0 y donde u < 0.
    Los logaritmos de `conservar` (los que ya traía el integrando, como ln x en ∫ x ln x dx) se dejan igual."""
    return F.replace(sp.log, lambda a: sp.log(a) if (_siempre_positivo(a) or isinstance(a, sp.Abs)
                                                     or sp.log(a) in conservar) else sp.log(sp.Abs(a)))


def sin_desfase(F):
    """sen(x + π/4) -> (sen x + cos x)/√2: simplify a veces junta un seno y un coseno en uno solo con desfase."""
    def es_desfase(e):
        return isinstance(e, (sp.sin, sp.cos)) and isinstance(e.args[0], sp.Add) and e.args[0].has(sp.pi)
    if not any(es_desfase(a) for a in F.atoms(sp.Function)):
        return F
    return sp.factor_terms(sp.expand(F.replace(es_desfase, sp.expand_trig)))


def sin_abs_en_log(f):
    """ln|u| -> ln u, solo para integrar (la derivada de las dos es u'/u); el valor absoluto se repone al final."""
    return f.replace(lambda e: isinstance(e, sp.log) and e.args[0].has(sp.Abs),
                     lambda e: sp.log(e.args[0].replace(sp.Abs, lambda a: a)))


def misma_funcion(g, h):
    """True si g y h coinciden: primero con sympy; si no decide, en los puntos donde las dos son reales."""
    try:
        if sp.simplify(g - h) == 0:
            return True
    except Exception:
        pass
    comparados = 0
    for c in _PUNTOS:
        try:
            if abs(valor_real(g, c) - valor_real(h, c)) > 1e-9 * max(1, abs(valor_real(h, c))):
                return False
            comparados += 1
        except (ValueError, TypeError, ZeroDivisionError, OverflowError):
            continue
    return comparados >= 3


def _a_logaritmos(F):
    """asinh, acosh y atanh se escriben con logaritmos (son elementales)."""
    return F.rewrite(sp.log) if F.has(sp.asinh, sp.acosh, sp.atanh) else F


def _rama_real(F, f):
    """De una antiderivada por casos (Piecewise) o con constantes imaginarias, una fórmula real cuya derivada es f."""
    candidatos = [e for e, _ in F.args] if isinstance(F, sp.Piecewise) else [F]
    for e in candidatos:
        for cand in (e, sp.re(e)):
            cand = _a_logaritmos(cand)
            if cand.has(sp.I, sp.re, sp.im, sp.Piecewise):
                continue
            if misma_funcion(derivada_legible(cand), f):
                return cand
    return None


def antiderivada(f):
    """Una antiderivada real de f (con C = 0), simplificada y con ln|..| donde corresponde."""
    g = sin_abs_en_log(f)
    F = con_limite(lambda: sp.integrate(g, x), 20)
    if F.has(sp.Integral):
        raise ValueError("sympy no encontró una antiderivada para esta función.")
    if F.has(sp.RootSum):
        raise ValueError("sympy expresó la antiderivada con raíces de un polinomio; este caso está fuera del curso.")
    F = sp.piecewise_fold(F)
    if F.has(sp.Piecewise, sp.I):
        F = _rama_real(F, g)
        if F is None:
            raise ValueError("sympy dio una antiderivada por casos que no pude reducir a una sola fórmula real.")
    F = _a_logaritmos(F)
    raras = {type(a).__name__ for a in F.atoms(sp.Function) if not isinstance(a, _PERMITIDAS)}
    if raras:
        raise ValueError(f"la antiderivada usa funciones especiales ({', '.join(sorted(raras))}): no es elemental. "
                         "Una integral definida de esta función se calcula con sumas numéricas (Unidad 4).")
    if F.is_polynomial(x):
        F = sp.expand(F)
    else:
        F = con_valor_absoluto(sin_desfase(_a_logaritmos(sp.trigsimp(sp.simplify(F)))), g.atoms(sp.log))
    if not misma_funcion(derivada_legible(F), g):
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


# Qué regla eligió sympy: nombres en español de las reglas de integral_steps
_NOMBRES = {
    "URule": "sustitución", "PartsRule": "integración por partes", "CyclicPartsRule": "partes (cíclica)",
    "RewriteRule": "reescritura del integrando (identidad, división o fracciones parciales)",
    "TrigSubstitutionRule": "sustitución trigonométrica",
    "SqrtQuadraticRule": "raíz de una cuadrática (sustitución trigonométrica)",
    "SqrtQuadraticDenomRule": "raíz de una cuadrática (sustitución trigonométrica)",
    "ReciprocalSqrtQuadraticRule": "raíz de una cuadrática (sustitución trigonométrica)",
    "CompleteSquareRule": "completar el cuadrado", "ArctanRule": "tabla (arco tangente)",
    "ArcsinRule": "tabla (arco seno)",
    "DontKnowRule": "integral_steps no reconoce la regla (sympy.integrate puede resolverla de otra forma)",
}
_TABLA = {"PowerRule", "ExpRule", "SinRule", "CosRule", "ReciprocalRule", "ConstantRule", "Sec2Rule"}
_ENVOLTURAS = {"AlternativeRule", "ConstantTimesRule"}


def _desenvuelve(r):
    while type(r).__name__ in _ENVOLTURAS:
        r = r.alternatives[0] if type(r).__name__ == "AlternativeRule" else r.substep
    return r


def regla_sympy(f):
    """Nombre en español de la técnica principal que integral_steps usa para f (la primera que no es envoltura)."""
    try:
        r = _desenvuelve(con_limite(lambda: _mi.integral_steps(f, x), 8))
    except ValueError:
        return "integral_steps tardó demasiado; se omite"
    nombre = type(r).__name__
    if nombre == "AddRule":
        partes = []
        for s_ in r.substeps:
            n_ = type(_desenvuelve(s_)).__name__
            partes.append("tabla de la Unidad 3" if n_ in _TABLA else _NOMBRES.get(n_, n_))
        return "suma término por término: " + ", ".join(dict.fromkeys(partes))
    if nombre in _TABLA:
        return "tabla de la Unidad 3"
    return _NOMBRES.get(nombre, nombre)


def _constante(dif):
    """Valor de una diferencia que debería ser constante: exacto si se puede, o el primero real de _PUNTOS."""
    dif = sp.simplify(dif)
    if not dif.has(x):
        v = dif
    else:
        v = None
        for c in _PUNTOS:
            try:
                v = valor_real(dif, c); break
            except (ValueError, TypeError, ZeroDivisionError, OverflowError):
                continue
        if v is None:
            return None
    vf = float(sp.N(v))
    if abs(vf) < 1e-10:
        return sp.Integer(0)
    return v if not isinstance(v, float) else sp.nsimplify(round(vf, 10), [sp.pi, sp.E], tolerance=1e-9)


def compara_con_sympy(F, f):
    """(es_correcta, constante): F' = f, y en cuánto difiere F de la antiderivada de sympy (None si no se sabe)."""
    if not misma_funcion(derivada_legible(F), sin_abs_en_log(f)):
        return False, None
    try:
        G = antiderivada(f)
    except ValueError:
        return True, None
    return True, _constante(F - G)


def mostrar_antiderivada(F):
    """F + C, con la constante siempre a la vista."""
    return f"{F} + C"


def texto_constante(cte):
    if cte is None:
        return ""
    return "coincide con la antiderivada de sympy" if cte == 0 else f"difiere de la antiderivada de sympy en la constante {cte}"
