# ID: ALG-U3-NB00
# Notebook: alg/u3_sistemas.ipynb · utilidades compartidas
# Repositorio: alg/u3_sistemas/00_utilidades.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# # Álg. U3 · Sistemas lineales y rango
# **Libro de regularización en matemáticas para ingeniería · ULSA Bajío**
#
# Este notebook acompaña la unidad 3 del bloque de álgebra lineal. Tiene una sección por subtema, en el mismo orden que el libro:
#
# 1. Eliminación de Gauss y Gauss-Jordan
# 2. Regla de Cramer
# 3. Ecuaciones vectorial y matricial; rango
# 4. Clasificación y sistemas homogéneos
# 5. Interpretación geométrica en 2D y 3D
#
# **Cómo usarlo.** Ejecuta las celdas en orden (Entorno de ejecución → Ejecutar todas). Cada sección tiene una celda de teoría, código que la recorre paso a paso, una calculadora, sus casos de prueba y una celda **Contrasta**: antes de ejecutarla, resuelve a mano el caso que se indica y escribe tu resultado. La referencia de la calculadora aparece solo después de que escribes el tuyo.
#
# **Etiqueta del repositorio.** `ALG-U3` (bloque, unidad). Las celdas de código de cada sección vienen de `alg/u3_sistemas/NN_*.py`; los fragmentos que el libro imprime, de `alg/u3_sistemas/fragmentos/`.
#
# **Convenciones.**
#
# - Punto decimal. Un sistema se escribe como su **matriz aumentada** $[A\,|\,\mathbf b]$, **por renglones**: cada renglón es una ecuación, los renglones se separan con `;` o con un salto de línea, y las entradas con comas o espacios. **La última columna es $\mathbf b$.** Por ejemplo, `1, 1, 6; 2, -1, 3` significa $x+y=6$, $2x-y=3$. También se acepta la forma de Python `[[1, 1, 6], [2, -1, 3]]` y una barra antes de $\mathbf b$: `1, 1 | 6; 2, -1 | 3`.
# - Puedes escribir fracciones como `2/3`; los decimales se convierten a fracciones exactas (`0.1` se lee como $\tfrac{1}{10}$), así que las cuentas internas no acumulan error de redondeo. Se acepta el signo menos tipográfico (−).
# - Las incógnitas se llaman $x, y, z$ cuando son tres o menos, y $x_1,\dots,x_6$ cuando son más. Las calculadoras aceptan hasta 6 ecuaciones y 6 incógnitas: una matriz aumentada de $6\times 7$ a lo más.
# - Los renglones se nombran $R_1, R_2,\dots$ y las columnas se cuentan desde 1. Las operaciones de renglón se escriben `R2 ← R2 - 2·R1` (al renglón 2 le restas 2 veces el renglón 1), `R2 ↔ R3` (intercambias los renglones 2 y 3) y `R1 ← (1/2)·R1` (multiplicas el renglón 1 por 1/2).
# - En lo que muestran las calculadoras, el signo menos es siempre `-`, también en las ecuaciones y en las operaciones de renglón.
# - Los decimales que se muestran se *redondean* (no se truncan) a las cifras significativas que elijas.

# %% [markdown]
# ## 0. Versiones e imports

# %%
# Versiones mínimas probadas. En Colab ya vienen instaladas; esta celda solo lo asegura.
%pip install -q "numpy>=1.26" "sympy>=1.12" "matplotlib>=3.8" "ipywidgets>=8"

import sys, math, re, io, textwrap, contextlib
from decimal import Decimal, ROUND_HALF_UP, localcontext
from fractions import Fraction
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registra la proyección "3d")
import sympy as sp
import ipywidgets as widgets
from IPython.display import display, Math

print("Python", sys.version.split()[0], "| numpy", np.__version__, "| sympy", sp.__version__,
      "| matplotlib", matplotlib.__version__, "| ipywidgets", widgets.__version__)

# %%
# Utilidades compartidas por todo el notebook: leer matrices exactas, redondear y mostrar
_NUMERO = re.compile(r"[+-]?(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?")
_AYUDA_NUMERO = "usa enteros, fracciones como 2/3 o decimales con punto"


_LIMITE = Fraction(10) ** 15                  # entradas entre 1e-15 y 1e15 en valor absoluto (o 0)


def _en_rango(q, texto):
    if q != 0 and not Fraction(1) / _LIMITE <= abs(q) <= _LIMITE:
        raise ValueError(f"«{texto}» está fuera del intervalo que maneja esta calculadora: "
                         "usa valores entre 1e-15 y 1e15 en valor absoluto (o 0).")
    return sp.Rational(q.numerator, q.denominator)


def numero(valor, nombre="el número"):
    """Un número como racional exacto de sympy: '0.1' -> 1/10, '2/3' -> 2/3, 5 -> 5, 1.2 -> 6/5."""
    if valor is None:
        raise ValueError(f"{nombre} está vacío: escribe un número.")
    if isinstance(valor, bool):
        raise ValueError(f"«{valor}» no es un número: {_AYUDA_NUMERO}.")
    if isinstance(valor, (int, np.integer)):
        return _en_rango(Fraction(int(valor)), valor)
    if isinstance(valor, (float, np.floating)) or (isinstance(valor, sp.Basic) and valor.is_Float):
        v = float(valor)
        if not math.isfinite(v):
            raise ValueError(f"{nombre} no es un número finito.")
        return _en_rango(Fraction(repr(v)), repr(v))      # repr da el decimal más corto: 0.1 -> '0.1' -> 1/10 exacto
    if isinstance(valor, sp.Basic):
        if valor.is_Rational:
            return valor                             # ya es exacto (por ejemplo, un resultado intermedio): no se acota
        raise ValueError(f"«{valor}» no es un número racional: {_AYUDA_NUMERO}.")
    texto = re.sub(r"\s+", "", str(valor)).replace("−", "-")      # acepta el signo menos tipográfico
    partes = texto.split("/")
    formas = [_NUMERO.fullmatch(p) for p in partes]
    if not texto or len(partes) > 2 or not all(formas):
        raise ValueError(f"«{str(valor).strip()}» no es un número: {_AYUDA_NUMERO}." if texto else
                         f"{nombre} está vacío: escribe un número.")
    if any(f.group(2) and abs(int(f.group(2)[1:])) > 30 for f in formas):   # 1e999999 no se calcula: se rechaza antes
        raise ValueError(f"«{texto}» está fuera del intervalo que maneja esta calculadora: "
                         "usa valores entre 1e-15 y 1e15 en valor absoluto (o 0).")
    num, *den = (Fraction(p) for p in partes)
    if den and den[0] == 0:
        raise ValueError(f"«{texto}» divide entre cero.")
    return _en_rango(num / den[0] if den else num, texto)


def _entradas(k):
    return f"{k} entrada" if k == 1 else f"{k} entradas"


def _desde_renglones(renglones, nombre):
    """Lista de renglones (cada uno, lista de entradas) -> Matrix de racionales, con mensajes claros."""
    if not renglones or not renglones[0]:
        raise ValueError(f"escribe {nombre} por renglones, por ejemplo: 1, 2; 3, 4")
    n = len(renglones[0])
    salida = []
    for i, renglon in enumerate(renglones, 1):
        if len(renglon) != n:
            raise ValueError(f"en {nombre}, el renglón {i} tiene {_entradas(len(renglon))} y el 1 tiene {n}: "
                             "todos los renglones deben tener el mismo número de entradas")
        fila = []
        for j, v in enumerate(renglon, 1):
            try:
                fila.append(numero(v))
            except ValueError as err:
                raise ValueError(f"en {nombre}, renglón {i}, columna {j}: {err}") from None
        salida.append(fila)
    return sp.Matrix(salida)


def matriz(texto, nombre="la matriz"):
    """Lee una matriz escrita por renglones ('1, 2; 3, 4' o un renglón por línea) como Matrix de racionales exactos."""
    if not isinstance(texto, str):
        return como_matriz(texto, nombre)
    t = texto.replace("−", "-")
    t = re.sub(r"[\]\)]\s*,?\s*[\[\(]", ";", t)            # [[1, 2], [3, 4]] o (1, 2), (3, 4): cada grupo es un renglón
    t = re.sub(r"[\[\]\(\)]", " ", t)
    renglones = [r.strip() for r in re.split(r"[;\n]", t)]
    renglones = [r for r in renglones if r]
    entradas = []
    for i, r in enumerate(renglones, 1):
        partes = re.split(r"\s*,\s*|\s+", r)
        if "" in partes:
            raise ValueError(f"en {nombre}, el renglón {i} tiene una entrada vacía: "
                             "revisa que no haya dos comas seguidas ni una coma al principio o al final")
        entradas.append(partes)
    try:
        return _desde_renglones(entradas, nombre)
    except ValueError as err:
        if "todos los renglones deben tener" in str(err):
            raise ValueError(f"{err}{pista_coma(texto)}") from None
        raise


def pista_coma(texto):
    """Si el texto tiene algo como 1,5 (coma entre dígitos, sin espacio), la pista de usar punto decimal."""
    return "; si escribiste decimales, usa punto decimal: 1.5, no 1,5" if isinstance(texto, str) and re.search(r"\d,\d", texto) else ""


def como_matriz(A, nombre="la matriz"):
    """Acepta texto, Matrix de sympy, lista de listas o arreglo de numpy y devuelve una Matrix de racionales exactos."""
    if isinstance(A, str):
        return matriz(A, nombre)
    if isinstance(A, sp.MatrixBase):
        renglones = A.tolist()
    elif isinstance(A, np.ndarray):
        if A.ndim > 2:
            raise ValueError(f"{nombre} debe ser una tabla de números (renglones y columnas)")
        renglones = np.atleast_2d(A).tolist()
    elif isinstance(A, (list, tuple)):
        anidada = [isinstance(f, (list, tuple, np.ndarray, sp.MatrixBase)) for f in A]
        if A and all(anidada):
            renglones = [list(f) for f in A]
        elif not any(anidada):
            renglones = [list(A)]                  # una lista sola es un renglón
        else:
            raise ValueError(f"en {nombre}, mezclaste números sueltos con renglones")
    else:
        renglones = [[A]]
    return _desde_renglones(renglones, nombre)


def vector(texto, nombre="el vector"):
    """Lee un vector ('1, 2, 3', '1; 2; 3', (1, 2, 3) o [1, 2, 3]) como columna de racionales exactos."""
    try:
        v = como_matriz(texto, nombre)
    except ValueError as err:
        if str(err).startswith("escribe "):          # vacío: el ejemplo debe ser de vector, no de matriz
            raise ValueError(f"escribe {nombre} como lista de números, por ejemplo: 1, 2, 3") from None
        raise
    if v.rows == 1 or v.cols == 1:
        return v.reshape(v.rows * v.cols, 1)
    raise ValueError(f"{nombre} debe ser una lista de números, por ejemplo: 1, 2, 3 (escribiste una matriz {dims(v)})")


def dims(M):
    """Tamaño 'm×n' de una matriz."""
    return f"{M.rows}×{M.cols}"


def cifras(valor, n=4):
    """Redondea (no trunca) a n cifras significativas con la regla escolar (5 sube) y conserva los ceros finales.
    Un racional exacto se redondea desde su valor exacto (p/q con 50 cifras), no desde su aproximación en punto
    flotante: 1249999999999999999/10^19 con 2 cifras da 0.12, no 0.13."""
    with localcontext() as ctx:
        ctx.prec = 50
        if isinstance(valor, (sp.Rational, Fraction, int, np.integer)) and not isinstance(valor, bool):
            p, q = (int(valor.p), int(valor.q)) if isinstance(valor, sp.Rational) else (
                (valor.numerator, valor.denominator) if isinstance(valor, Fraction) else (int(valor), 1))
            d = Decimal(p) / Decimal(q)
        else:
            v = float(valor)
            if math.isnan(v):
                return "no definido"
            if math.isinf(v):
                return "+∞" if v > 0 else "-∞"
            d = Decimal(repr(v))
        if d == 0:
            return "0"
        q = d.quantize(Decimal(1).scaleb(d.adjusted() - n + 1), rounding=ROUND_HALF_UP)
        if q.adjusted() != d.adjusted():      # 9.996 con 3 cifras pasa a 10.0: se reajusta la posición
            q = d.quantize(Decimal(1).scaleb(q.adjusted() - n + 1), rounding=ROUND_HALF_UP)
        e = q.adjusted()
        if e > 15 or e < -4:                  # muy grandes o muy chicos: notación científica, 1.234×10^-7
            return f"{format(q.scaleb(-e), 'f')}×10^{e}"
        if q.adjusted() >= n:                 # enteros grandes: sin notación científica
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


def _mensaje(fn):
    """Texto del ValueError que lanza fn() ('' si no lanza): para probar que el mensaje explica el error."""
    try:
        fn()
    except ValueError as err:
        return str(err)
    return ""


# --- formato de números y matrices ---------------------------------------------------------------
_SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")


def _sub(*indices):
    """Subíndices: _sub(1, 2) -> '₁₂'; con índices de dos cifras los separa con coma."""
    if all(k < 10 for k in indices):
        return "".join(str(k) for k in indices).translate(_SUB)
    return "₍" + ",".join(str(k) for k in indices).translate(_SUB) + "₎"


def _paren(v):
    """Un número listo para ir en un producto: negativos y fracciones entre paréntesis."""
    v = sp.sympify(v)
    return str(v) if (v >= 0 and v.is_integer) else f"({v})"


def _en_linea(M):
    """Matriz en una línea con la misma notación de entrada: [2, 3; 5, -2]."""
    return "[" + "; ".join(", ".join(str(e) for e in M.row(i)) for i in range(M.rows)) + "]"


def _tupla(v, n=None):
    """Vector como (a, b, c); con n, además sus decimales redondeados a n cifras."""
    exacto = "(" + ", ".join(str(e) for e in v) + ")"
    if n is None or all(e.is_integer for e in v):
        return exacto
    return exacto + " ≈ (" + ", ".join(cifras(e, n) for e in v) + ")"


def _aprox(v, n):
    """Decimales redondeados de un vector: '(0.5000, 1.333)'; con una sola entrada, sin paréntesis."""
    textos = [cifras(e, n) for e in v]
    return textos[0] if len(textos) == 1 else "(" + ", ".join(textos) + ")"


def _plano(M, fmt=str, corte=None):
    """Texto plano de una matriz con columnas alineadas; corte = columna tras la que se dibuja una barra."""
    celdas = [[fmt(M[i, j]) for j in range(M.cols)] for i in range(M.rows)]
    anchos = [max(len(celdas[i][j]) for i in range(M.rows)) for j in range(M.cols)]
    lineas = []
    for fila in celdas:
        partes = [c.rjust(a) for c, a in zip(fila, anchos)]
        if corte:
            partes = partes[:corte] + ["|"] + partes[corte:]
        lineas.append("[ " + "  ".join(partes) + " ]")
    return "\n".join(lineas)


_TEX = [("⁻¹", "^{-1}"), ("ᵀ", "^{T}"), ("·", r"\cdot "), ("×", r"\times "), ("≈", r"\approx "),
        ("−", "-"), ("det ", r"\det "), ("|", r"\,|\,")]
_SUB_TEX = re.compile("[₀₁₂₃₄₅₆₇₈₉]+")


_PALABRAS = re.compile(r"(?<![\\A-Za-zÁÉÍÓÚáéíóúñ])([A-Za-zÁÉÍÓÚáéíóúñ]{3,}(?:\s+[A-Za-zÁÉÍÓÚáéíóúñ]{2,})*)")


def _a_tex(nombre):
    """Nombre de una cantidad en LaTeX: 'det A₁' -> '\\det A_{1}'; las palabras ('forma escalonada de') van en \\text."""
    for a, b in _TEX:
        nombre = nombre.replace(a, b)
    nombre = _SUB_TEX.sub(lambda m: "_{" + m.group(0).translate(str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")) + "}", nombre)
    return _PALABRAS.sub(lambda m: r"\text{" + m.group(1) + "} ", nombre)


def _tex_decimales(M, n):
    filas = r" \\ ".join(" & ".join(cifras(M[i, j], n) for j in range(M.cols)) for i in range(M.rows))
    return r"\left[\begin{matrix}" + filas + r"\end{matrix}\right]"


class _Formula(Math):
    """Math de IPython con una versión en texto plano (la que se ve donde no se dibuja LaTeX)."""
    def __init__(self, tex, texto):
        super().__init__(tex)
        self._texto = texto

    def __repr__(self):
        return self._texto


class _Bloque:
    """Salida en Markdown (texto en recuadro y fórmulas en LaTeX) con su versión en texto plano."""
    def __init__(self, partes):
        self._partes = partes

    @staticmethod
    def _envuelve(texto, ancho=100):
        """Corta los renglones largos (en un recuadro de código no se ajustan solos); conserva la sangría."""
        lineas = []
        for linea in texto.rstrip("\n").split("\n"):
            sangria = " " * (len(linea) - len(linea.lstrip(" ")))
            lineas += textwrap.wrap(linea, ancho, subsequent_indent=sangria + "  ") if len(linea) > ancho else [linea]
        return "\n".join(lineas)

    def _repr_markdown_(self):
        return "\n\n".join(f"$${tex}$$" if tipo == "tex" else "```text\n" + self._envuelve(texto) + "\n```"
                            for tipo, tex, texto in self._partes)

    def __repr__(self):
        return "\n".join(texto.rstrip("\n") for _, _, texto in self._partes)


_BLOQUES = []      # bloques abiertos: lo que se muestra dentro de uno se junta y se entrega de una vez


@contextlib.contextmanager
def bloque():
    """Junta lo que se imprime y las fórmulas que se muestran dentro del bloque y lo entrega como UNA sola salida.
    Las calculadoras lo usan para que su resultado aparezca completo de una vez, en el orden en que se produjo."""
    partes, buf = [], io.StringIO()

    def vacia():
        if buf.getvalue():
            partes.append(("texto", "", buf.getvalue()))
            buf.seek(0); buf.truncate()

    _BLOQUES.append((partes, vacia))
    try:
        with contextlib.redirect_stdout(buf):
            yield
    finally:
        vacia()
        _BLOQUES.pop()
        if partes:
            _emite(_Bloque(partes)) if _BLOQUES else display(_Bloque(partes))


def _emite(formula):
    """Muestra una fórmula, o la agrega al bloque abierto si lo hay."""
    if not _BLOQUES:
        display(formula)
        return
    partes, vacia = _BLOQUES[-1]
    vacia()
    if isinstance(formula, _Bloque):
        partes.extend(formula._partes)
    else:
        partes.append(("tex", formula.data, formula._texto))


def muestra(nombre, M, n=None):
    """Muestra 'nombre = M' con LaTeX (con su tamaño si es matriz). Con n agrega los decimales redondeados a n cifras."""
    if isinstance(M, sp.MatrixBase):
        tex = f"{_a_tex(nombre)} = {sp.latex(M)}"
        texto = f"{nombre}  ({dims(M)}) =\n{_plano(M)}"
        if n is not None and not all(e.is_integer for e in M):
            tex += r" \approx " + _tex_decimales(M, n)
            texto += f"\n≈\n{_plano(M, lambda e: cifras(e, n))}"
        tex += rf"\qquad ({M.rows}\times {M.cols})"
    else:
        M = sp.sympify(M)
        tex, texto = f"{_a_tex(nombre)} = {sp.latex(M)}", f"{nombre} = {M}"
        if n is not None and M.is_number and not M.is_integer:
            tex += rf" \approx {cifras(M, n)}"
            texto += f" ≈ {cifras(M, n)}"
    _emite(_Formula(tex, texto))


def compara(mia, ref, pista=""):
    """'coinciden' o en qué entradas difiere tu matriz de la referencia, más una pista."""
    if mia.shape != ref.shape:
        return f"NO coinciden: tu matriz es {dims(mia)} y la de referencia es {dims(ref)}. {pista}".strip()
    malas = [f"({i + 1}, {j + 1})" for i in range(ref.rows) for j in range(ref.cols) if mia[i, j] != ref[i, j]]
    if not malas:
        return "coinciden"
    cuales = "la entrada" if len(malas) == 1 else "las entradas"
    return f"NO coinciden en {cuales} {', '.join(malas)} (renglón, columna). {pista}".strip()


# %%
# Utilidades nuevas de esta unidad: sistemas como matriz aumentada, eliminación exacta y solución general
_MAX_ECUACIONES, _MAX_INCOGNITAS = 6, 6       # matriz aumentada de 6×7 a lo más


def _ecuaciones(k):
    return f"{k} ecuación" if k == 1 else f"{k} ecuaciones"


def _incognitas(k):
    return f"{k} incógnita" if k == 1 else f"{k} incógnitas"


def _revisa_tamano(m, n, nombre, cosa="incógnitas"):
    """Rechaza más de 6 renglones o más de 6 incógnitas (columnas de A)."""
    limite = "una matriz aumentada de 6×7 a lo más" if cosa == "incógnitas" else "matrices de hasta 6×6"
    if m > _MAX_ECUACIONES:
        raise ValueError(f"{nombre} tiene {m} renglones; esta calculadora acepta hasta {_MAX_ECUACIONES} ({limite})")
    if n > _MAX_INCOGNITAS:
        if cosa == "incógnitas":
            raise ValueError(f"{nombre} tiene {n + 1} columnas, es decir, {n} incógnitas más la columna b; esta calculadora "
                             f"acepta hasta {_MAX_INCOGNITAS} incógnitas ({limite})")
        raise ValueError(f"{nombre} tiene {n} columnas; esta calculadora acepta hasta {_MAX_INCOGNITAS} ({limite})")


def aumentada(texto, nombre="la matriz aumentada"):
    """Lee un sistema escrito como su matriz aumentada [A | b] por renglones; la última columna es b.
    '1, 1, 6; 2, -1, 3' es x + y = 6, 2x - y = 3. Acepta una barra antes de b: '1, 1 | 6; 2, -1 | 3'."""
    if isinstance(texto, str):
        texto = texto.replace("|", ",")
    try:
        M = como_matriz(texto, nombre)
    except ValueError as err:
        if str(err).startswith("escribe "):          # vacío: el ejemplo debe ser de un sistema
            raise ValueError(f"escribe {nombre} por renglones, por ejemplo: 1, 1, 6; 2, -1, 3 "
                             "(cada renglón es una ecuación y la última columna es b)") from None
        raise
    if M.cols < 2:
        raise ValueError(f"{nombre} necesita al menos dos columnas: los coeficientes y, al final, la columna b "
                         "(por ejemplo, «2, 6» es la ecuación 2x = 6)")
    try:
        _revisa_tamano(M.rows, M.cols - 1, nombre)
    except ValueError as err:
        raise ValueError(f"{err}{pista_coma(texto)}") from None
    return M


def matriz_limitada(texto, nombre="A"):
    """Una matriz (no aumentada) de hasta 6×6, con los mismos mensajes que matriz()."""
    M = como_matriz(texto, nombre)
    try:
        _revisa_tamano(M.rows, M.cols, nombre, cosa="columnas")
    except ValueError as err:
        raise ValueError(f"{err}{pista_coma(texto)}") from None
    return M


def de_columnas(*columnas):
    """Matriz cuyas columnas son los vectores dados: de_columnas('1, 0, 1', '0, 1, 1')."""
    vs = [vector(c, f"la columna {j}") for j, c in enumerate(columnas, 1)]
    if not vs:
        raise ValueError("escribe al menos una columna")
    if len({v.rows for v in vs}) > 1:
        raise ValueError("todas las columnas deben tener el mismo número de entradas")
    return sp.Matrix.hstack(*vs)


def nombres_incognitas(n):
    """x, y, z con tres incógnitas o menos; x₁, ..., xₙ con más."""
    return ["x", "y", "z"][:n] if n <= 3 else [f"x{_sub(k)}" for k in range(1, n + 1)]


def parametros(k):
    """Nombres de k parámetros, en el orden de las variables libres: s, t, u, v, w (con más de cinco: s₁, s₂, ...)."""
    if k == 0:
        return []
    if k <= 5:
        return list(sp.symbols("s t u v w")[:k])
    return list(sp.symbols(" ".join(f"s{_sub(i)}" for i in range(1, k + 1))))


def _factor(a):
    """Coeficiente delante de un renglón: 1 -> '', 2 -> '2·', 1/2 -> '(1/2)·', -1 -> '(-1)·'."""
    if a == 1:
        return ""
    if a.is_integer and a > 0:
        return f"{a}·"
    return f"({a})·"


def _op_suma(i, m, k):
    """Texto de R_i ← R_i - m·R_k con el signo que toca: m = 2 -> 'R2 ← R2 - 2·R1'; m = -3 -> 'R2 ← R2 + 3·R1'."""
    signo, a = ("-", m) if m > 0 else ("+", -m)
    return f"R{i} ← R{i} {signo} {_factor(a)}R{k}"


def _elimina(M, n_a, con_b=True):
    """Eliminación de Gauss y después de Jordan, exacta, sobre una copia de M (Matrix de racionales).

    Los pivotes se buscan en las primeras n_a columnas (las de A); se intercambian renglones solo si el pivote es 0
    y abajo hay una entrada distinta de cero (como a mano). Con con_b, al final se busca un pivote en la columna n_a
    (la de b): si aparece, hay un renglón 0 = c y la forma escalonada lo es también de [A | b].
    La fase de Jordan deja cada pivote en 1 y ceros arriba de él, también el de la columna de b si lo hay: en un sistema
    incompatible, el renglón 0 = c queda como 0 = 1 y la columna de b se limpia hacia arriba (forma reducida de verdad).
    'imposible' guarda el renglón y el valor c de la forma escalonada, antes de normalizar.
    Devuelve un diccionario: E (escalonada), R (reducida), pivotes (columnas de A con pivote, desde 0), libres
    (columnas de A sin pivote), r (rango de A), r_aum (rango de [A | b]), imposible ((renglón, c) o None) y pasos,
    una lista de (tipo, texto, matriz después de la operación o None) con tipo 'fase', 'nota' u 'op'."""
    E = M.copy()
    m = E.rows
    pasos = [("fase", "Gauss: ceros debajo de cada pivote (forma escalonada)", None)]
    pivotes, libres, f, imposible = [], [], 0, None
    for c in range(n_a + (1 if con_b and M.cols > n_a else 0)):
        es_b = c == n_a
        k = next((i for i in range(f, m) if E[i, c] != 0), None)
        if k is None:
            if not es_b:
                libres.append(c)
                motivo = "ya no quedan renglones abajo" if f >= m else f"de R{f + 1} hacia abajo todo es cero"
                pasos.append(("nota", f"la columna {c + 1} no tiene pivote: {motivo}", None))
            continue
        if k != f:
            pasos.append(("nota", f"el pivote de la columna {c + 1} es cero; se intercambian R{f + 1} y R{k + 1}", None))
            E.row_swap(f, k)
            pasos.append(("op", f"R{f + 1} ↔ R{k + 1}", E.copy()))
        for i in range(f + 1, m):
            if E[i, c] != 0:
                mult = E[i, c] / E[f, c]
                E[i, :] = E[i, :] - mult * E[f, :]
                pasos.append(("op", _op_suma(i + 1, mult, f + 1), E.copy()))
        if es_b:
            imposible = (f + 1, E[f, c])
            pasos.append(("nota", f"R{f + 1} quedó como 0 = {E[f, c]}: ningún valor de las incógnitas lo cumple", None))
        else:
            pivotes.append(c)
        f += 1
    if len(pasos) == 1:
        pasos.append(("nota", "la matriz ya estaba escalonada: no hace falta ninguna operación", None))
    R = E.copy()
    pasos.append(("fase", "Jordan: cada pivote en 1 y ceros también arriba (forma escalonada reducida)", None))
    hechas = len(pasos)
    todos = list(enumerate(pivotes)) + ([(len(pivotes), n_a)] if imposible is not None else [])
    for f, c in reversed(todos):
        if R[f, c] != 1:
            k = 1 / R[f, c]
            R[f, :] = k * R[f, :]
            pasos.append(("op", f"R{f + 1} ← {_factor(k)}R{f + 1}", R.copy()))
        for i in range(f):
            if R[i, c] != 0:
                mult = R[i, c]
                R[i, :] = R[i, :] - mult * R[f, :]
                pasos.append(("op", _op_suma(i + 1, mult, f + 1), R.copy()))
    if len(pasos) == hechas:
        pasos.append(("nota", "la forma escalonada ya era reducida: no hace falta ninguna operación", None))
    return {"E": E, "R": R, "pivotes": pivotes, "libres": libres, "r": len(pivotes),
            "r_aum": len(pivotes) + (imposible is not None), "imposible": imposible, "pasos": pasos}


def operaciones(el):
    """Solo los textos de las operaciones de renglón de un resultado de _elimina."""
    return [texto for tipo, texto, _ in el["pasos"] if tipo == "op"]


def _imprime_pasos(pasos, corte=None, sangria="  "):
    for tipo, texto, M in pasos:
        if tipo == "fase":
            print(f"— {texto} —")
        elif tipo == "nota":
            print(f"{sangria}Nota: {texto}.")
        else:
            print(f"{sangria}{texto}")
            print("\n".join(sangria + "  " + linea for linea in _plano(M, corte=corte).split("\n")))


def solucion_general(R, pivotes, n):
    """Solución de un sistema compatible a partir de su forma reducida R (con la columna b al final).
    Devuelve (x, particular, direcciones, parámetros): cada variable libre es un parámetro (s, t, u, ... en el orden
    de las libres), x = particular + s·d₁ + t·d₂ + ..., y la particular es la que tiene las libres en 0."""
    libres = [c for c in range(n) if c not in pivotes]
    ps = parametros(len(libres))
    particular = sp.zeros(n, 1)
    for f, c in enumerate(pivotes):
        particular[c] = R[f, n]
    direcciones = []
    for c_libre in libres:
        d = sp.zeros(n, 1)
        d[c_libre] = 1
        for f, c in enumerate(pivotes):
            d[c] = -R[f, c_libre]
        direcciones.append(d)
    x = particular + sum((p * d for p, d in zip(ps, direcciones)), sp.zeros(n, 1))
    return x, particular, direcciones, ps


def _coeficiente(a, letra, primero):
    """Un término 'a·letra' con el signo que toca: (-2, 's', False) -> '- 2s'; (1/4, 's', True) -> '(1/4)s'."""
    mag = abs(a)
    coef = "" if (mag == 1 and letra) else (f"{mag}" if mag.is_integer else f"({mag})")
    if primero:
        return f"{'-' if a < 0 else ''}{coef}{letra}"
    return f"{'-' if a < 0 else '+'} {coef}{letra}"


def _lineal(e, ps):
    """Expresión lineal en los parámetros, con la constante primero: '4 - 2s', '(1/4)s', '7 - 2s - 3t'."""
    e = sp.expand(sp.sympify(e))
    constante = e.subs({p: 0 for p in ps})
    partes = [str(constante)] if constante != 0 else []
    for p in ps:
        a = e.coeff(p)
        if a != 0:
            partes.append(_coeficiente(a, str(p), not partes))
    return " ".join(partes) if partes else "0"


def _ecuacion(fila, nombres):
    """Renglón [a₁, ..., aₙ, b] como ecuación: '2x - y + z = 3'; si todos los coeficientes son 0, '0 = b'."""
    terminos = []
    for a, v in zip(list(fila)[:-1], nombres):
        if a != 0:
            terminos.append(_coeficiente(a, v, not terminos))
    return f"{' '.join(terminos) if terminos else '0'} = {list(fila)[-1]}"


def _imprime_sistema(Ab, sangria="  "):
    m, n = Ab.rows, Ab.cols - 1
    print(f"Sistema de {_ecuaciones(m)} con {_incognitas(n)}:")
    for i in range(m):
        print(f"{sangria}{_ecuacion(Ab.row(i), nombres_incognitas(n))}")


def _tex_aumentada(M, corte):
    columnas = "c" * corte + "|" + "c" * (M.cols - corte)
    filas = r" \\ ".join(" & ".join(sp.latex(M[i, j]) for j in range(M.cols)) for i in range(M.rows))
    return r"\left[\begin{array}{" + columnas + "}" + filas + r"\end{array}\right]"


def muestra_aumentada(nombre, M, corte):
    """Muestra una matriz aumentada con la barra después de la columna 'corte'."""
    _emite(_Formula(f"{_a_tex(nombre)} = {_tex_aumentada(M, corte)}", f"{nombre} =\n{_plano(M, corte=corte)}"))


def _lista(cosas):
    """Números en una frase: [1, 2, 3] -> '1, 2 y 3'. (Para nombres de incógnitas usa _coma: 'x y z' confundiría.)"""
    cosas = [str(c) for c in cosas]
    return cosas[0] if len(cosas) == 1 else ", ".join(cosas[:-1]) + " y " + cosas[-1]


def _coma(cosas):
    """Nombres de incógnitas o parámetros separados con comas: ['x', 'y'] -> 'x, y'."""
    return ", ".join(str(c) for c in cosas)


def _flotante(M):
    """Matrix de sympy -> arreglo de numpy en punto flotante."""
    return np.array(M.tolist(), dtype=float)


# Gráficas: un color y un estilo de línea por ecuación, en orden fijo, para que se lean también en blanco y negro
_COLORES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
_ESTILOS = ["-", "--", ":", "-.", (0, (5, 1, 1, 1, 1, 1)), (0, (1, 4))]          # continuo, discontinuo, punteado...


_LLAMADAS_SOLVE = [0]      # cuántas veces se llamó a numpy.linalg.solve (las pruebas revisan que no se llame de más)


_AVISO_SINGULAR = ("numpy.linalg.solve no pudo resolverlo en punto flotante (la matriz es singular para la máquina), "
                   "aunque la solución exacta existe: el sistema está muy mal condicionado")


def _numpy_solve(A, b):
    """(x, None) con numpy.linalg.solve en punto flotante, solo para comparar con la solución exacta; (None, aviso)
    si numpy no puede: una matriz invertible en fracciones puede ser singular para la máquina."""
    _LLAMADAS_SOLVE[0] += 1
    try:
        with np.errstate(all="ignore"):
            x = np.linalg.solve(_flotante(A), _flotante(b).ravel())
    except np.linalg.LinAlgError:
        return None, _AVISO_SINGULAR
    if not np.all(np.isfinite(x)):
        return None, _AVISO_SINGULAR
    return x, None


def normaliza_texto(texto):
    """Minúsculas, sin acentos ni signos de puntuación y con un solo espacio entre palabras: '¿Se cortan?' -> 'se cortan'."""
    import unicodedata
    t = unicodedata.normalize("NFD", str(texto).lower())
    t = "".join(c for c in t if unicodedata.category(c) != "Mn")
    return " ".join(re.sub(r"[^a-z0-9ñ ]+", " ", t).split())


def compara_vector(mio, ref, nombres, pista=""):
    """'coinciden' o qué incógnitas difieren de la referencia, más una pista."""
    if mio.rows != ref.rows:
        return (f"NO coinciden: escribiste {_entradas(mio.rows).replace('entrada', 'número')} y el sistema tiene "
                f"{_incognitas(ref.rows)}. {pista}").strip()
    malas = [nombres[k] for k in range(ref.rows) if mio[k] != ref[k]]
    if not malas:
        return "coinciden"
    return f"NO coinciden en {_coma(malas)}. {pista}".strip()


def evalua(pruebas):
    """Ejecuta una vez cada caso (nombre, función que devuelve True o False). Devuelve [(nombre, ok)]; un caso
    que lanza una excepción cuenta como falla y su nombre dice por qué, sin interrumpir a los demás."""
    resultados = []
    for nombre, prueba in pruebas:
        try:
            ok = bool(prueba())
        except Exception as err:
            ok, nombre = False, f"{nombre}  (error: {type(err).__name__}: {err})"
        resultados.append((nombre, ok))
    return resultados


def corre_pruebas(pruebas, seccion):
    """Escribe OK o FALLA en cada caso y el resumen de la sección."""
    resultados = evalua(pruebas)
    for nombre, ok in resultados:
        print("OK  " if ok else "FALLA", nombre)
    print(f"Sección {seccion}: {sum(ok for _, ok in resultados)}/{len(pruebas)} casos pasan")
