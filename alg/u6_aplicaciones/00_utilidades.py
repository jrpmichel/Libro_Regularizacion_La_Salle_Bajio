# ID: ALG-U6-NB00
# Notebook: alg/u6_aplicaciones.ipynb · utilidades compartidas
# Repositorio: alg/u6_aplicaciones/00_utilidades.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# # Álg. U6 · Aplicaciones a ingeniería
# **Libro de regularización en matemáticas para ingeniería · ULSA Bajío**
#
# Este notebook acompaña la unidad 6 del bloque de álgebra lineal. Cada sección resuelve con vectores, matrices y sistemas un problema que vas a encontrar en las materias de tu carrera, en el mismo orden que el libro:
#
# 1. Equilibrio de fuerzas en 3D, momento y trabajo
# 2. Análisis de mallas con leyes de Kirchhoff
# 3. Mezclas y asignación de recursos
# 4. Armadura por el método de nodos
# 5. Ajuste de datos por mínimos cuadrados
# 6. Transformaciones de coordenadas
#
# **Cómo usarlo.** Ejecuta las celdas en orden (Entorno de ejecución → Ejecutar todas). Cada sección tiene una celda de teoría, código que la recorre paso a paso, una calculadora, sus casos de prueba y una celda **Contrasta**: antes de ejecutarla, resuelve a mano el caso que se indica y escribe tu resultado. La referencia de la calculadora aparece solo después de que escribes el tuyo.
#
# **Etiqueta del repositorio.** `ALG-U6` (bloque, unidad). Las celdas de código de cada sección vienen de `alg/u6_aplicaciones/NN_*.py`; los fragmentos que el libro imprime, de `alg/u6_aplicaciones/fragmentos/`.
#
# **Convenciones.**
#
# - **Las listas se escriben un elemento por renglón** (o con los elementos separados por `;`), y los números de cada elemento se separan con comas o espacios. Cada calculadora dice qué lleva cada renglón: un anclaje `x, y, z` (sección 1); una rama `R, malla, malla` y una fuente `V, malla` (sección 2); un ingrediente con sus porcentajes (sección 3); un nodo `nombre, x, y`, una barra `nombre, nombre`, un apoyo `nombre, x` o `nombre, y` y una carga `nombre, Fx, Fy` (sección 4).
# - Un vector suelto (una carga, un punto, una traslación) se escribe como lista de números: `0, 0, -90`. También se acepta `(0, 0, -90)`.
# - **Los datos de la sección 5** se escriben como listas de números separados por espacios o comas: `0 2 3.6` o `0, 2, 3.6` (hasta 500).
# - Punto decimal. Puedes escribir fracciones como `2/3` (también `2 / 3`); los decimales se convierten a fracciones exactas (`0.3` se lee como $\tfrac{3}{10}$), así que las tensiones, las corrientes, las cantidades de una mezcla y las fuerzas en las barras se calculan sin error de redondeo. Las raíces que salen de las longitudes (por ejemplo $\sqrt2$ en un vector unitario) se dejan exactas y se muestran también en decimales. Se aceptan los signos menos tipográficos (− y –) y un espacio después del signo (`- 1`).
# - **Unidades.** Cada calculadora dice en qué unidades van sus datos (o tiene un campo para escribirlas) y las escribe en sus resultados. Un ángulo pide **siempre** su unidad, `"grados"` o `"radianes"`: no hay unidad por omisión.
# - Mallas, nodos de un circuito, cables, ingredientes y datos se cuentan desde 1 (los pasos de tiempo de la sección 5, desde 0: $x[0],x[1],\dots$). En lo que muestran las calculadoras, el signo menos es siempre `-` y las raíces se escriben con `√`. Los decimales que se muestran se *redondean* (no se truncan) a las cifras significativas que elijas; si una fracción exacta resulta demasiado larga, la calculadora muestra solo sus decimales y lo avisa.
# - Si usas las funciones desde código, les puedes pasar texto (como en las calculadoras) o listas de Python: `equilibrio_particula([(2, 2, 1), (-2, 1, 2), (1, -2, 2)], (0, 0, -90))`.

# %% [markdown]
# ## 0. Versiones e imports

# %%
# Versiones mínimas probadas. En Colab ya vienen instaladas; esta celda solo lo asegura.
%pip install -q "numpy>=1.26" "sympy>=1.12" "matplotlib>=3.8" "ipywidgets>=8"

import sys, math, re, io, textwrap, contextlib, time
from decimal import Decimal, ROUND_HALF_UP, localcontext
from fractions import Fraction
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401  (registra la proyección "3d")
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import sympy as sp
import ipywidgets as widgets
from IPython.display import display, Math

print("Python", sys.version.split()[0], "| numpy", np.__version__, "| sympy", sp.__version__,
      "| matplotlib", matplotlib.__version__, "| ipywidgets", widgets.__version__)

# %%
# Utilidades compartidas por todo el notebook (las mismas de las unidades 4 y 5): leer números y matrices
# exactos, redondear y mostrar
matplotlib.rcParams["axes.unicode_minus"] = False     # signo menos ASCII también en los ejes de las gráficas
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
            raise ValueError(f"{nombre} no es finito (escribiste {v}).")
        return _en_rango(Fraction(repr(v)), repr(v))      # repr da el decimal más corto: 0.1 -> '0.1' -> 1/10 exacto
    if isinstance(valor, sp.Basic):
        if valor.is_Rational:
            return valor                             # ya es exacto (por ejemplo, un resultado intermedio): no se acota
        raise ValueError(f"«{valor}» no es un número racional: {_AYUDA_NUMERO}.")
    texto = re.sub(r"\s+", "", str(valor)).replace("−", "-").replace("–", "-")      # acepta los signos menos tipográficos
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


def _componentes(k):
    return f"{k} componente" if k == 1 else f"{k} componentes"


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


def _normaliza_entrada(texto):
    """Signos menos tipográficos (− y –) como '-', llaves como paréntesis y sin espacios alrededor de '/' ni entre un
    signo y su número: '2 / 3' -> '2/3', '- 1' -> '-1', '{(1, 1), (0, 1)}' -> '((1, 1), (0, 1))'."""
    t = texto.replace("−", "-").replace("–", "-").replace("{", "(").replace("}", ")")
    t = re.sub(r"\s*/\s*", "/", t)
    return re.sub(r"(^|[\s,;(\[])([+-])\s+(?=[\d.])", r"\1\2", t)


def _trozos(texto):
    """Texto -> lista de renglones, cada uno lista de entradas (texto). Acepta ';' o saltos de línea entre renglones,
    comas o espacios entre entradas, y la forma de Python [[1, 2], [3, 4]], (1, 2), (3, 4) o {(1, 2), (3, 4)}.
    Devuelve también el número del primer renglón con una entrada vacía (o None)."""
    t = _normaliza_entrada(texto)
    t = re.sub(r"[\]\)]\s*,?\s*[\[\(]", ";", t)            # [[1, 2], [3, 4]] o (1, 2), (3, 4): cada grupo es un renglón
    t = re.sub(r"[\[\]\(\)]", " ", t)
    renglones = [r.strip() for r in re.split(r"[;\n]", t)]
    renglones = [r for r in renglones if r]
    entradas, vacia = [], None
    for i, r in enumerate(renglones, 1):
        partes = re.split(r"\s*,\s*|\s+", r)
        if "" in partes and vacia is None:
            vacia = i
        entradas.append(partes)
    return entradas, vacia


def matriz(texto, nombre="la matriz"):
    """Lee una matriz escrita por renglones ('1, 2; 3, 4' o un renglón por línea) como Matrix de racionales exactos."""
    if not isinstance(texto, str):
        return como_matriz(texto, nombre)
    entradas, vacia = _trozos(texto)
    if vacia is not None:
        raise ValueError(f"en {nombre}, el renglón {vacia} tiene una entrada vacía: "
                         "revisa que no haya dos comas seguidas ni una coma al principio o al final")
    try:
        return _desde_renglones(entradas, nombre)
    except ValueError as err:
        if "todos los renglones deben tener" in str(err):
            raise ValueError(f"{err}{pista_coma(texto)}") from None
        raise


def pista_coma(texto):
    """La pista de usar punto decimal, solo si la lectura es de veras ambigua: hay comas entre dígitos sin espacio
    (1,5) y además números separados solo por espacios (1,5 2,5), que es como se escriben decimales con coma."""
    if isinstance(texto, str) and re.search(r"\d,\d", texto) and re.search(r"\d[ \t]+[-+.]?\d", texto):
        return "; si escribiste decimales, usa punto decimal: 1.5, no 1,5"
    return ""


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
    """Lee un vector ('1, 2, 3', '1; 2; 3', (1, 2, 3) o [1, 2, 3]) como columna de racionales exactos.
    Los mensajes de error dicen qué componente está mal."""
    ejemplo = f"escribe {nombre} como lista de números, por ejemplo: 5, 4, 1"
    if texto is None:
        raise ValueError(ejemplo)
    if isinstance(texto, str):
        filas, vacia = _trozos(texto)
        if vacia is not None:
            raise ValueError(f"en {nombre} hay una componente vacía: revisa que no haya dos comas seguidas "
                             "ni una coma al principio o al final")
        if not filas:
            raise ValueError(ejemplo)
        if len(filas) == 1:
            entradas = filas[0]
        elif all(len(f) == 1 for f in filas):
            entradas = [f[0] for f in filas]
        else:
            raise ValueError(f"{nombre} debe ser una lista de números, por ejemplo: 5, 4, 1 (escribiste varios renglones)")
    elif isinstance(texto, (list, tuple, np.ndarray, sp.MatrixBase)):
        plano = ((isinstance(texto, np.ndarray) and texto.ndim == 1) or
                 (isinstance(texto, (list, tuple)) and
                  not any(isinstance(e, (list, tuple, np.ndarray, sp.MatrixBase)) for e in texto)))
        if not plano:                                   # Matrix, arreglo 2D o lista de listas
            return _vector_de_matriz(texto, nombre, ejemplo)
        entradas = list(texto)
        if not entradas:
            raise ValueError(ejemplo)
    else:
        entradas = [texto]                              # un número solo
    valores = []
    for j, e in enumerate(entradas, 1):
        try:
            valores.append(numero(e))
        except ValueError as err:
            raise ValueError(f"en {nombre}, componente {j}: {err}") from None
    return sp.Matrix(valores)


def _vector_de_matriz(texto, nombre, ejemplo):
    """Un vector dado como Matrix, arreglo o lista anidada: debe tener un solo renglón o una sola columna."""
    try:
        v = como_matriz(texto, nombre)
    except ValueError as err:
        if str(err).startswith("escribe "):
            raise ValueError(ejemplo) from None
        raise
    if v.rows == 1 or v.cols == 1:
        return v.reshape(v.rows * v.cols, 1)
    raise ValueError(f"{nombre} debe ser una lista de números, por ejemplo: 5, 4, 1 (escribiste una matriz {dims(v)})")


def dims(M):
    """Tamaño 'm×n' de una matriz."""
    return f"{M.rows}×{M.cols}"


def cifras(valor, n=4):
    """Redondea (no trunca) a n cifras significativas con la regla escolar (5 sube) y conserva los ceros finales.
    Un racional exacto se redondea desde su valor exacto (p/q con 50 cifras), no desde su aproximación en punto
    flotante: 1249999999999999999/10^19 con 2 cifras da 0.12, no 0.13. Una expresión exacta con raíces, senos o
    cosenos se evalúa con 50 cifras antes de redondear."""
    with localcontext() as ctx:
        ctx.prec = 50
        if isinstance(valor, sp.Basic) and not valor.is_Rational and valor.is_number:
            valor = sp.N(valor, 50)
            if valor.is_real is False:
                return "no definido"
            d = Decimal(str(valor)) if valor.is_Float else Decimal(repr(float(valor)))
        elif isinstance(valor, (sp.Rational, Fraction, int, np.integer)) and not isinstance(valor, bool):
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
        if q.adjusted() >= n:                 # enteros grandes: exactos tal cual; redondeados, con ×10^
            return f"{int(q):d}" if q == d else f"{format(q.scaleb(-e), 'f')}×10^{e}"
        return format(q, "f")


def decimales_fijos(valor, d=2):
    """Redondea (5 sube) a d decimales fijos, desde el valor exacto: sqrt(3) con 2 -> '1.73'; 1 -> '1.00'."""
    with localcontext() as ctx:
        ctx.prec = 50
        x = Decimal(str(sp.N(sp.sympify(valor), 40)))
        return format(x.quantize(Decimal(1).scaleb(-d), rounding=ROUND_HALF_UP), "f")


def _rechaza(fn):
    """True si fn() lanza ValueError (para probar que una función rechaza una entrada)."""
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
_SUP = str.maketrans("0123456789", "⁰¹²³⁴⁵⁶⁷⁸⁹")


def _sub(*indices):
    """Subíndices: _sub(1, 2) -> '₁₂'; con índices de dos cifras los separa con coma."""
    if all(k < 10 for k in indices):
        return "".join(str(k) for k in indices).translate(_SUB)
    return "₍" + ",".join(str(k) for k in indices).translate(_SUB) + "₎"


def _R(n):
    """Nombre del espacio en texto plano: _R(3) -> 'R³'."""
    return "R" + str(n).translate(_SUP)


def _txt(e):
    """Una cantidad exacta como texto legible: sqrt(3)/2 -> √3/2, 2*sin(30) -> 2·sen(30), pi/6 -> π/6."""
    t = str(e)
    t = re.sub(r"sqrt\((\d+)\)", r"√\1", t).replace("sqrt(", "√(")
    return t.replace("**", "^").replace("*", "·").replace("sin(", "sen(").replace("pi", "π")


def _latex(x):
    """LaTeX de sympy con la convención del libro: sen en lugar de sin."""
    return sp.latex(x).replace(r"\sin", r"\operatorname{sen}")


def _exacta_legible(M, largo=28, total=400):
    """True si la forma exacta de M se lee bien: cada entrada con pocos caracteres y no demasiadas en total.
    Si no, las calculadoras muestran solo los decimales y lo avisan."""
    textos = [_txt(e) for e in M]
    return all(len(t) <= largo for t in textos) and sum(len(t) for t in textos) <= total


_LARGA = "la forma exacta es muy larga: se muestra redondeada"


def _igual(nombre, x, n=4):
    """'nombre = x' con x exacto (y sus decimales si no es entero); si la forma exacta es muy larga,
    'nombre ≈ decimales (aviso)'."""
    x = sp.sympify(x)
    if len(_txt(x)) > 28:
        return f"{nombre} ≈ {cifras(x, n)} ({_LARGA})"
    return f"{nombre} = {_txt(x)}" + ("" if x.is_integer else f" ≈ {cifras(x, n)}")


def _en_linea(M):
    """Matriz en una línea con la misma notación de entrada: [2, 3; 5, -2]."""
    return "[" + "; ".join(", ".join(_txt(e) for e in M.row(i)) for i in range(M.rows)) + "]"


def _tupla(v, n=None):
    """Vector como (a, b, c); con n, además sus decimales redondeados a n cifras. Si la forma exacta es muy larga,
    solo los decimales (a 4 cifras si no se da n), con el aviso."""
    if not _exacta_legible(v):
        return "≈ (" + ", ".join(cifras(e, n or 4) for e in v) + f") ({_LARGA})"
    exacto = "(" + ", ".join(_txt(e) for e in v) + ")"
    if n is None or all(e.is_integer for e in v):
        return exacto
    return exacto + " ≈ (" + ", ".join(cifras(e, n) for e in v) + ")"


def _aprox(v, n):
    """Decimales redondeados de un vector: '(0.5000, 1.333)'; con una sola entrada, sin paréntesis."""
    textos = [cifras(e, n) for e in v]
    return textos[0] if len(textos) == 1 else "(" + ", ".join(textos) + ")"


def _plano(M, fmt=None, corte=None):
    """Texto plano de una matriz con columnas alineadas; corte = columna tras la que se dibuja una barra."""
    fmt = fmt or _txt
    celdas = [[fmt(M[i, j]) for j in range(M.cols)] for i in range(M.rows)]
    anchos = [max(len(celdas[i][j]) for i in range(M.rows)) for j in range(M.cols)]
    lineas = []
    for fila in celdas:
        partes = [c.rjust(a) for c, a in zip(fila, anchos)]
        if corte:
            partes = partes[:corte] + ["|"] + partes[corte:]
        lineas.append("[ " + "  ".join(partes) + " ]")
    return "\n".join(lineas)


_TEX = [("⁻¹", "^{-1}"), ("ᵀ", "^{T}"), ("²", "^{2}"), ("·", r"\cdot "), ("×", r"\times "), ("≈", r"\approx "),
        ("←", r"\leftarrow "), ("−", "-"), ("|", r"\,|\,")]
_SUB_TEX = re.compile("[₀₁₂₃₄₅₆₇₈₉]+")
_PROSA = re.compile(r"[A-Za-zÁÉÍÓÚáéíóúñü]{3,}")


def _a_tex(nombre):
    """Nombre de una cantidad en LaTeX. Si lleva palabras ('forma escalonada reducida'), va completo en \\text{...},
    con sus espacios. Si no, se traducen los símbolos (P_C⁻¹·P_B -> P_C^{-1}\\cdot P_B) y los subíndices (A2 -> A_{2},
    v1 -> \\mathbf v_{1}). Las calculadoras pasan el LaTeX a mano cuando el nombre mezcla palabras y símbolos."""
    if _PROSA.search(re.sub(r"\bdet\b", "", nombre)):
        return r"\text{" + nombre.replace("{", r"\{").replace("}", r"\}") + "}"
    t = nombre
    for a, b in _TEX:
        t = t.replace(a, b)
    t = re.sub(r"\bdet\s+", r"\\det ", t)
    t = _SUB_TEX.sub(lambda m: "_{" + m.group(0).translate(str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")) + "}", t)
    return re.sub(r"(?<![A-Za-z0-9\\])([A-Za-z])(\d+)(?![\d/.])",
                  lambda m: (r"\mathbf " if m.group(1) in "vbeu" else "") + m.group(1) + "_{" + m.group(2) + "}", t)


def _tex_cifra(t):
    """Un número ya redondeado, en LaTeX: '1.234×10^-7' -> '1.234\\times 10^{-7}'."""
    return re.sub(r"×10\^(-?\d+)", r"\\times 10^{\1}", t)


def _tex_decimales(M, n):
    filas = r" \\ ".join(" & ".join(_tex_cifra(cifras(M[i, j], n)) for j in range(M.cols)) for i in range(M.rows))
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
            partes.append(("texto", "", buf.getvalue().replace("= ≈", "≈")))     # un valor que salió redondeado
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


def muestra(nombre, M, n=None, tex=None):
    """Muestra 'nombre = M' con LaTeX (con su tamaño si es matriz). Con n agrega los decimales redondeados a n cifras.
    tex es el nombre en LaTeX (si no se da, se arma con _a_tex). Si la forma exacta es muy larga, solo los decimales."""
    tex_nombre = tex if tex is not None else _a_tex(nombre)
    if isinstance(M, sp.MatrixBase):
        if _exacta_legible(M):
            tex = f"{tex_nombre} = {_latex(M)}"
            texto = f"{nombre}  ({dims(M)}) =\n{_plano(M)}"
            if n is not None and not all(e.is_integer for e in M):
                tex += r" \approx " + _tex_decimales(M, n)
                texto += f"\n≈\n{_plano(M, lambda e: cifras(e, n))}"
        else:
            tex = f"{tex_nombre} \\approx {_tex_decimales(M, n or 4)}" + r"\quad\text{(" + _LARGA + ")}"
            texto = f"{nombre}  ({dims(M)}) ≈\n{_plano(M, lambda e: cifras(e, n or 4))}\n({_LARGA})"
        tex += rf"\qquad ({M.rows}\times {M.cols})"
    else:
        M = sp.sympify(M)
        tex, texto = f"{tex_nombre} = {_latex(M)}", f"{nombre} = {_txt(M)}"
        if n is not None and M.is_number and not M.is_integer:
            tex += rf" \approx {_tex_cifra(cifras(M, n))}"
            texto += f" ≈ {cifras(M, n)}"
    _emite(_Formula(tex, texto))


# %%
# Utilidades de esta unidad: tablas con un elemento por renglón, unidades, eliminación exacta, gráficas, pruebas
# y lectura de listas de datos
def _es_cero(e):
    """True si e vale 0. Un racional se compara directo; un número exacto (con raíces, senos o cosenos) se evalúa con
    50 cifras, sin simplificarlo (simplificar radicales anidados puede tardar minutos); una expresión con símbolos se
    simplifica."""
    e = sp.sympify(e)
    if e.is_Rational:
        return e == 0
    if e.is_number:                                   # en decimales con 50 cifras: rápido y suficiente
        return bool(abs(sp.N(e, 50)) < sp.Float("1e-40", 50))
    return sp.simplify(e) == 0


def _simplifica_uno(e):
    """Una entrada simplificada. Un número que vale 0, 1 o -1 (comprobado con 50 cifras) se escribe así; uno con senos o
    cosenos, o muy largo (más de 40 operaciones), se deja como está: simplify lo reconstruiría y sympy podría
    desarrollar el coseno en radicales anidados, lo que puede tardar minutos."""
    if e.is_Rational:
        return e
    if e.is_number:
        for v in (0, 1, -1):
            if _es_cero(e - v):
                return sp.Integer(v)
        if e.has(sp.sin, sp.cos) or sp.count_ops(e) > 40:
            return e
    return sp.simplify(e)


def _simplifica(M):
    """Simplifica las entradas que no son racionales (las racionales se quedan como están), con _simplifica_uno."""
    return M.applyfunc(_simplifica_uno)


def _factor(a):
    """Coeficiente delante de un renglón: 1 -> '', 2 -> '2·', 1/2 -> '(1/2)·', -1 -> '(-1)·'."""
    if a == 1:
        return ""
    if a.is_integer and a > 0:
        return f"{a}·"
    return f"({a})·"


def _op_suma(i, m, k):
    """Texto de R_i ← R_i - m·R_k con el signo que toca: m = 2 -> 'R2 ← R2 - 2·R1'; m = -3 -> 'R2 ← R2 + 3·R1'."""
    if m.is_number and m.is_real and m < 0:
        return f"R{i} ← R{i} + {_factor(-m)}R{k}"
    return f"R{i} ← R{i} - {_factor(m)}R{k}"


def _elimina(M, n_a, con_b=True):
    """Eliminación de Gauss y después de Jordan, exacta, sobre una copia de M (la misma de las unidades 4 y 5).

    Los pivotes se buscan en las primeras n_a columnas; se intercambian renglones solo si el pivote es 0 y abajo hay
    una entrada distinta de cero (como a mano). Con con_b, al final se busca un pivote en la columna n_a (la de b):
    si aparece, hay un renglón 0 = c. Devuelve un diccionario: E (escalonada), R (reducida), pivotes (columnas con
    pivote, desde 0), libres, r (rango), r_aum (rango con la columna b), imposible ((renglón, c) o None) y pasos."""
    E = M.copy()
    racional = all(e.is_Rational for e in E)
    limpia = (lambda X: X) if racional else _simplifica
    m = E.rows
    pasos = [("fase", "Gauss: ceros debajo de cada pivote (forma escalonada)", None)]
    pivotes, libres, f, imposible = [], [], 0, None
    for c in range(n_a + (1 if con_b and M.cols > n_a else 0)):
        es_b = c == n_a
        k = next((i for i in range(f, m) if not _es_cero(E[i, c])), None)
        if k is None:
            if not es_b:
                libres.append(c)
            continue
        if k != f:
            E.row_swap(f, k)
            pasos.append(("op", f"R{f + 1} ↔ R{k + 1}", E.copy()))
        for i in range(f + 1, m):
            if not _es_cero(E[i, c]):
                mult = limpia(sp.Matrix([E[i, c] / E[f, c]]))[0]
                E[i, :] = limpia(E[i, :] - mult * E[f, :])
                E[i, c] = 0
                pasos.append(("op", _op_suma(i + 1, mult, f + 1), E.copy()))
        if es_b:
            imposible = (f + 1, E[f, c])
        else:
            pivotes.append(c)
        f += 1
    R = E.copy()
    pasos.append(("fase", "Jordan: cada pivote en 1 y ceros también arriba (forma escalonada reducida)", None))
    todos = list(enumerate(pivotes)) + ([(len(pivotes), n_a)] if imposible is not None else [])
    for f, c in reversed(todos):
        if not _es_cero(R[f, c] - 1):
            k = limpia(sp.Matrix([1 / R[f, c]]))[0]
            R[f, :] = limpia(k * R[f, :])
            R[f, c] = 1
            pasos.append(("op", f"R{f + 1} ← {_factor(k)}R{f + 1}", R.copy()))
        for i in range(f):
            if not _es_cero(R[i, c]):
                mult = R[i, c]
                R[i, :] = limpia(R[i, :] - mult * R[f, :])
                R[i, c] = 0
                pasos.append(("op", _op_suma(i + 1, mult, f + 1), R.copy()))
    return {"E": E, "R": R, "pivotes": pivotes, "libres": libres, "r": len(pivotes),
            "r_aum": len(pivotes) + (imposible is not None), "imposible": imposible, "pasos": pasos}


def _direcciones(R, pivotes, n):
    """Base de Nul(A) leída de la forma reducida R de A (sus primeras n columnas): un vector por variable libre,
    con esa libre en 1 y las demás libres en 0."""
    direcciones = []
    for c_libre in [c for c in range(n) if c not in pivotes]:
        d = sp.zeros(n, 1)
        d[c_libre] = 1
        for f, c in enumerate(pivotes):
            d[c] = -R[f, c_libre]
        direcciones.append(d)
    return direcciones


def base_nulo(M):
    """Base exacta de Nul(M) = {x : M x = 0} (lista vacía si solo contiene al 0)."""
    el = _elimina(M, M.cols, con_b=False)
    return _direcciones(el["R"], el["pivotes"], M.cols)


def _lista(cosas):
    """Números en una frase: [1, 2, 3] -> '1, 2 y 3'."""
    cosas = [str(c) for c in cosas]
    return cosas[0] if len(cosas) == 1 else ", ".join(cosas[:-1]) + " y " + cosas[-1]


def _coma(cosas):
    """Nombres separados con comas: ['v1', 'v2'] -> 'v1, v2'."""
    return ", ".join(str(c) for c in cosas)


def _flotante(M):
    """Matrix de sympy -> arreglo de numpy en punto flotante."""
    return np.array(sp.Matrix(M).evalf().tolist(), dtype=float)


def normaliza_texto(texto):
    """Minúsculas, sin acentos ni signos de puntuación y con un solo espacio entre palabras: '¿Dirección?' -> 'direccion'."""
    import unicodedata
    t = unicodedata.normalize("NFD", str(texto).lower())
    t = "".join(c for c in t if unicodedata.category(c) != "Mn")
    return " ".join(re.sub(r"[^a-z0-9ñ ]+", " ", t).split())


# Gráficas: un color, un estilo de línea y un marcador por elemento, en orden fijo, para que se lean también en
# blanco y negro
_COLORES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
_ESTILOS = ["-", "--", ":", "-.", (0, (5, 1, 1, 1.5, 1, 1.5)), (0, (1, 4))]          # continuo, discontinuo, punteado...
_MARCAS = ["o", "s", "^", "D", "v", "P"]


def _flecha(ax, origen, punta, color, estilo, ancho, texto=None, desplazamiento=(6, 6)):
    """Flecha en el plano de origen a punta; con texto, su rótulo junto a la punta (sobre un fondo blanco)."""
    if np.hypot(*(np.asarray(punta, float) - np.asarray(origen, float))) > 1e-12:
        ax.annotate("", xy=punta, xytext=origen,
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=ancho, linestyle=estilo, shrinkA=0, shrinkB=0))
    if texto:
        ax.annotate(texto, xy=punta, xytext=desplazamiento, textcoords="offset points", fontsize=10,
                    bbox=dict(boxstyle="round,pad=0.15", facecolor="white", edgecolor="none", alpha=0.85))


def _captura(fn, *args, **kwargs):
    """Ejecuta fn sin mostrar nada (ni texto ni gráficas) y devuelve (texto, error): lo que imprimió o mostró y la
    excepción que lanzó (None si no lanzó). Sirve para las pruebas de humo de las calculadoras."""
    partes, buf = [], io.StringIO()

    def vacia():
        if buf.getvalue():
            partes.append(("texto", "", buf.getvalue()))
            buf.seek(0); buf.truncate()

    mostrar, error = plt.show, None
    _BLOQUES.append((partes, vacia))
    plt.show = lambda *a, **k: plt.close("all")
    try:
        with contextlib.redirect_stdout(buf):
            fn(*args, **kwargs)
    except Exception as err:
        error = err
    finally:
        vacia()
        _BLOQUES.pop()
        plt.show = mostrar
        plt.close("all")
    return "\n".join(texto for _, _, texto in partes), error


def _humo(calculadora, controles, **cambios):
    """Prueba de humo: la calculadora corre con los valores de sus controles (y los cambios dados), no lanza ninguna
    excepción, escribe algo y no dice «Revisa la entrada»."""
    valores = {k: w.value for k, w in controles.items()}
    valores.update(cambios)
    texto, error = _captura(calculadora, **valores)
    return error is None and texto.strip() != "" and "Revisa la entrada" not in texto


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


def _cerca(a, b, tol=1e-9):
    """True si |a - b| <= tol (para comparar resultados en decimales en las pruebas)."""
    return abs(float(a) - float(b)) <= tol


_COND_AVISO = 1000               # número de condición a partir del cual se avisa que la matriz es casi singular


def _condicion(M, por_renglones=False):
    """Número de condición (en decimales) de una matriz exacta: grande si es casi singular; inf si es singular. Con
    por_renglones=True, cada renglón se divide antes entre su norma, para que renglones de escalas muy distintas (un
    renglón de unos y otro de porcentajes chicos) no den un número grande sin que la matriz sea casi singular."""
    F = _flotante(M)
    if por_renglones:
        normas = np.linalg.norm(F, axis=1)
        normas[normas == 0] = 1.0
        F = F / normas[:, None]
    with np.errstate(all="ignore"):
        c = float(np.linalg.cond(F))
    return c if math.isfinite(c) else math.inf


class Diagnostico(ValueError):
    """Un hallazgo físico legítimo (cables en un mismo plano, Z singular, ingredientes iguales...), no un error de
    formato: las calculadoras lo muestran como «Diagnóstico: ...», y los errores de formato como «Revisa la entrada: ...»."""


def _reporta(err):
    """Lo que escribe una calculadora cuando no puede calcular: «Diagnóstico: ...» si es un hallazgo físico y
    «Revisa la entrada: ...» si es un error en lo que se escribió."""
    print("Diagnóstico:" if isinstance(err, Diagnostico) else "Revisa la entrada:", err)


def _punto(u, v):
    """Producto interior exacto de dos vectores columna: u·v = uᵀv."""
    return (u.T * v)[0]


def _es_nulo(v):
    """True si v es el vector cero."""
    return all(_es_cero(e) for e in v)


def _valor(x, n=4):
    """Un valor exacto como texto, con sus decimales si no es entero: '2', '1/3 ≈ 0.3333', '6·√2 ≈ 8.485'. Si la forma
    exacta es muy larga, solo los decimales: '≈ 1.234 (aviso)'. Un valor decimal (de numpy) se redondea: '≈ 0.8'."""
    if isinstance(x, (float, np.floating)):
        return f"≈ {cifras(x, n)}"
    x = sp.sympify(x)
    if len(_txt(x)) > 28:
        return f"≈ {cifras(x, n)} ({_LARGA})"
    return _txt(x) + ("" if x.is_integer else f" ≈ {cifras(x, n)}")


def _eq(nombre, texto):
    """'nombre = texto', o 'nombre ≈ ...' si el texto ya es aproximado (empieza con '≈')."""
    return f"{nombre} {texto}" if texto.startswith("≈") else f"{nombre} = {texto}"


def _decimal_exacto(e, n=15):
    """El texto decimal EXACTO de un racional con decimal finito de a lo más n cifras significativas ('0.3', '2.5',
    '-80', '0.0199'), o None si no lo tiene (1/3, o más de n cifras)."""
    if isinstance(e, (int, np.integer)) and not isinstance(e, bool):
        e = sp.Integer(int(e))
    if not isinstance(e, sp.Rational):
        return None
    if e == 0:
        return "0"
    q = int(e.q)
    for f in (2, 5):
        while q % f == 0:
            q //= f
    if q != 1:
        return None
    with localcontext() as ctx:
        ctx.prec = 60
        d = (Decimal(int(e.p)) / Decimal(int(e.q))).normalize()
    if len(d.as_tuple().digits) > n or not -6 <= d.adjusted() <= 15:
        return None
    return format(d, "f")


def _dec_o_cifras(x, n=4):
    """Para rótulos de gráficas: el decimal exacto si es corto ('30', '0.8', '-2.5'); si no, redondeado a n cifras."""
    if isinstance(x, sp.Basic) and x.is_Rational:
        corto = _decimal_exacto(x, n)
        if corto is not None:
            return corto
    return cifras(x, n)


def _tupla_dec(v, n=4):
    """Vector para rótulos de gráficas: '(2, 0.5, 0)', '(0.7071, 0.7071, 0)'."""
    return "(" + ", ".join(_dec_o_cifras(sp.sympify(e), n) for e in v) + ")"


def _ntxt(e, n=15, cif=4):
    """Un valor exacto como texto, en decimal si su decimal es finito y corto ('0.3', '2.5') y si no, con fracciones
    y raíces ('1/3', '6·√2'). Si la forma exacta es muy larga, sus decimales redondeados a cif cifras con el aviso."""
    corto = _decimal_exacto(e, n)
    if corto is not None:
        return corto
    t = _txt(sp.sympify(e))
    return t if len(t) <= 28 else f"≈ {cifras(e, cif)} ({_LARGA})"


def _magnitud(nombre, x, unidad, n=4):
    """Una cantidad con su unidad: 'T1 = 30 N', 'i1 = 2.5 A' (decimal exacto de a lo más n cifras), 'T1 = 6·√2 N ≈ 8.485 N'
    o 'AB = 20/3 kN ≈ 6.667 kN' (decimales redondeados a n cifras). Si la forma exacta es muy larga o el valor viene de
    numpy: 'R ≈ 20.17 Ω'."""
    u = f" {unidad}" if unidad else ""
    if isinstance(x, (float, np.floating)):
        return f"{nombre} ≈ {cifras(x, n)}{u}"
    x = sp.sympify(x)
    corto = _decimal_exacto(x, n)
    if corto is not None:
        return f"{nombre} = {corto}{u}"
    if len(_txt(x)) > 28:
        return f"{nombre} ≈ {cifras(x, n)}{u} ({_LARGA})"
    return f"{nombre} = {_txt(x)}{u} ≈ {cifras(x, n)}{u}"


def _combinacion(coefs, letra):
    """Lado izquierdo de una ecuación: '6·i1 - 4·i2', '0.46·x1 + 0.18·x2'. Solo los términos con coeficiente distinto
    de cero; el coeficiente se escribe siempre, también el 1."""
    partes = []
    for j, z in enumerate(coefs, 1):
        if z == 0:
            continue
        signo = ("-" if z < 0 else "") if not partes else ("- " if z < 0 else "+ ")
        partes.append(f"{signo}{_ntxt(abs(z))}·{letra}{j}")
    return " ".join(partes) if partes else "0"


def _cantidad(x, unidad, n=4):
    """Un valor con su unidad, sin nombre: '30 N', '6·√2 N ≈ 8.485 N', '≈ 20.17 Ω'."""
    t = _magnitud("", x, unidad, n)
    return t[3:] if t.startswith(" = ") else t[1:]


def _cuenta(k, singular, plural):
    """'1 dato', '5 datos': el número con el sustantivo en singular o en plural."""
    return f"{k} {singular if k == 1 else plural}"


def _aviso_coma(*textos):
    """Aviso (aunque no haya error) si la entrada tiene comas entre dígitos y números separados solo por espacios,
    que es como se escriben los decimales con coma: la calculadora leyó cada coma como separador."""
    if any(pista_coma(t) for t in textos):
        print("Aviso: hay comas entre dígitos (1,5) y números separados solo por espacios. Si escribiste decimales "
              "con coma, usa punto decimal (1.5, no 1,5): la calculadora leyó cada coma como separador.")


def nombre_unidad(texto, que, ejemplo):
    """La unidad que escribe el alumno (N, kN, kg, L, s...), sin espacios a los lados. Vacía o muy larga: ValueError."""
    t = str(texto if texto is not None else "").strip()
    if not t:
        raise ValueError(f"escribe la unidad de {que}, por ejemplo {ejemplo}")
    if len(t) > 12:
        raise ValueError(f"la unidad de {que} es muy larga (escribiste «{t}»): usa un símbolo corto, por ejemplo {ejemplo}")
    return t


def _sin_articulo(cosa):
    """'la rama' -> 'rama', 'el nodo' -> 'nodo'."""
    partes = cosa.split(" ", 1)
    return partes[1] if len(partes) == 2 and partes[0] in ("el", "la") else cosa


def _indefinido(cosa):
    """'la rama' -> 'una rama', 'el nodo' -> 'un nodo' (el artículo viene con el sustantivo: concordancia de género)."""
    return ("una " if cosa.startswith("la ") else "un ") + _sin_articulo(cosa)


def _pista_decimal(texto):
    """La pista del punto decimal si hay una coma entre dos dígitos (3,6), para listas donde la coma separa datos."""
    if isinstance(texto, str) and re.search(r"\d,\d", texto):
        return "; si escribiste decimales, usa punto decimal: 3.6, no 3,6"
    return ""


def _filas_crudas(texto, nombre, cosa, ejemplo):
    """Texto o lista de Python -> lista de renglones, cada uno lista de entradas (texto o número), sin convertir.
    cosa lleva su artículo ('la rama', 'el nodo'), para que los mensajes concuerden en género."""
    if texto is None:
        raise ValueError(f"escribe {nombre} con {_indefinido(cosa)} por renglón, por ejemplo: {ejemplo}")
    if isinstance(texto, str):
        filas, vacia = _trozos(texto)
        if vacia is not None:
            raise ValueError(f"en {nombre}, el renglón {vacia} tiene una entrada vacía: revisa que no haya dos comas "
                             "seguidas ni una coma al principio o al final")
    elif isinstance(texto, sp.MatrixBase):
        filas = texto.tolist()
    elif isinstance(texto, np.ndarray):
        filas = np.atleast_2d(texto).tolist()
    elif isinstance(texto, (list, tuple)):
        filas = []
        for i, f in enumerate(texto, 1):
            if isinstance(f, str):
                sub, vacia = _trozos(f)
                if vacia is not None or len(sub) != 1:
                    raise ValueError(f"en {nombre}, {cosa} {i} debe ir en un solo renglón, por ejemplo: {ejemplo}")
                filas.append(sub[0])
            elif isinstance(f, (list, tuple, np.ndarray, sp.MatrixBase)):
                filas.append(list(f))
            else:
                raise ValueError(f"en {nombre}, escribe cada {_sin_articulo(cosa)} como una lista, por ejemplo: {ejemplo}")
    else:
        raise ValueError(f"escribe {nombre} con {_indefinido(cosa)} por renglón, por ejemplo: {ejemplo}")
    if not filas or not any(filas):
        raise ValueError(f"escribe {nombre} con {_indefinido(cosa)} por renglón, por ejemplo: {ejemplo}")
    return filas


def tabla(texto, nombre, cosa, columnas, maximo, ejemplo):
    """Lee una lista de elementos (anclajes, ramas, fuentes, ingredientes...) escrita UN ELEMENTO POR RENGLÓN o con
    los elementos separados por ';', con `columnas` números cada uno (None: todos con los mismos que el primero).
    También acepta una lista de Python [(2, 1, 0), (4, 1, 2)] (las ramas de la sección 2). Devuelve una lista de
    renglones de racionales exactos. cosa lleva su artículo ('la rama', 'el anclaje'). Los mensajes de error dicen qué
    elemento y qué número están mal."""
    filas = _filas_crudas(texto, nombre, cosa, ejemplo)
    if len(filas) > maximo:
        raise ValueError(f"{nombre}: escribiste {len(filas)} renglones y esta calculadora acepta hasta {maximo} "
                         f"(recuerda: cada renglón es {_indefinido(cosa)})")
    k = columnas if columnas is not None else len(filas[0])
    salida = []
    for i, f in enumerate(filas, 1):
        if len(f) != k:
            raise ValueError(f"en {nombre}, {cosa} {i} tiene {_cuenta(len(f), 'número', 'números')} y debe tener "
                             f"{k}, por ejemplo: {ejemplo}{pista_coma(texto) or _pista_decimal(texto)}")
        renglon = []
        for j, e in enumerate(f, 1):
            try:
                renglon.append(numero(e))
            except ValueError as err:
                raise ValueError(f"en {nombre}, {_sin_articulo(cosa)} {i}, número {j}: {err}") from None
        salida.append(renglon)
    return salida


def _entero(valor, donde, minimo, maximo):
    """Un número entero entre minimo y maximo (un número de malla o de nodo): ValueError claro si no lo es."""
    if not valor.is_integer or not minimo <= valor <= maximo:
        raise ValueError(f"{donde} debe ser un entero de {minimo} a {maximo} (escribiste {_txt(valor)})")
    return int(valor)


def _tupla_u(v, unidad, n=None):
    """Vector con su unidad: '(0.3, 0.4, 0) m' (decimales exactos), '(-80, 60, 0) N·m'; si alguna componente no tiene
    decimal finito y corto, en forma exacta y, con n, también sus decimales redondeados: '(5/3, 1/3) A ≈ (1.667, 0.3333) A'."""
    u = f" {unidad}" if unidad else ""
    if not _exacta_legible(v):
        return "≈ (" + ", ".join(cifras(e, n or 4) for e in v) + f"){u} ({_LARGA})"
    cortos = [_decimal_exacto(e, n or 15) for e in v]
    if all(c is not None for c in cortos):
        return "(" + ", ".join(cortos) + ")" + u
    exacto = "(" + ", ".join(_txt(e) for e in v) + ")" + u
    if n is None:
        return exacto
    return exacto + " ≈ (" + ", ".join(cifras(e, n) for e in v) + ")" + u


_MAX_DATOS = 500                 # datos de la sección 5: hasta 500 por lista
_MAX_EXACTO = 12                 # con hasta 12 ecuaciones (o datos) el cálculo es exacto; con más, en decimales


def _texto_dato(e):
    """Un dato como se escribió: el texto tal cual; un número de Python o de sympy, en su forma corta."""
    if isinstance(e, str):
        return e
    if isinstance(e, (float, np.floating)):
        return repr(float(e))
    if isinstance(e, sp.Basic):
        return _txt(e)
    return str(e)


def datos(texto, nombre="los datos", con_texto=False):
    """Lee una lista de números (separados por espacios, comas, ';' o saltos de línea) como lista de racionales
    exactos: '0 2 3.6' -> [0, 2, 18/5]. Acepta también una lista o un arreglo de números. Hasta 500 datos; los
    mensajes de error dicen qué dato está mal. Con con_texto=True devuelve también los datos como se escribieron."""
    ejemplo = f"escribe {nombre} como lista de números separados por espacios o comas, por ejemplo: 0 1 2 3"
    if texto is None:
        raise ValueError(ejemplo)
    if isinstance(texto, str):
        filas, vacia = _trozos(texto)
        if vacia is not None:
            raise ValueError(f"en {nombre} hay un dato vacío: revisa que no haya dos comas seguidas ni una coma al "
                             "principio o al final")
        entradas = [e for fila in filas for e in fila]
    elif isinstance(texto, np.ndarray):
        entradas = list(texto.ravel())
    elif isinstance(texto, sp.MatrixBase):
        entradas = list(texto)
    elif isinstance(texto, (list, tuple)):
        if any(isinstance(e, (list, tuple, np.ndarray, sp.MatrixBase)) for e in texto):
            raise ValueError(f"{nombre} debe ser una lista de números, sin listas dentro")
        entradas = list(texto)
    else:
        entradas = [texto]                              # un número solo
    if not entradas:
        raise ValueError(ejemplo)
    if len(entradas) > _MAX_DATOS:
        raise ValueError(f"{nombre} son {len(entradas)} números; esta calculadora acepta hasta {_MAX_DATOS}")
    valores = []
    for j, e in enumerate(entradas, 1):
        try:
            valores.append(numero(e))
        except ValueError as err:
            raise ValueError(f"en {nombre}, dato {j}: {err}") from None
    return (valores, [_texto_dato(e) for e in entradas]) if con_texto else valores


def _num(c, n=4):
    """Un número para tablas y ecuaciones: si es un racional exacto con decimal finito de a lo más n cifras
    significativas, ese decimal tal cual (10.02, 1.5, -0.2, 3); si no, redondeado a n cifras con cifras()."""
    if isinstance(c, (int, np.integer)) and not isinstance(c, bool):
        c = sp.Integer(int(c))
    if isinstance(c, sp.Rational):
        q = int(c.q)
        for f in (2, 5):
            while q % f == 0:
                q //= f
        if q == 1:
            with localcontext() as ctx:
                ctx.prec = 60
                d = (Decimal(int(c.p)) / Decimal(int(c.q))).normalize()
            if d == 0:
                return "0"
            if len(d.as_tuple().digits) <= n and -4 <= d.adjusted() <= 15:
                return format(d, "f")
    return cifras(c, n)


def muestra_decimales(nombre, M, n=4, tex=None):
    """Muestra una matriz o un vector de punto flotante (de numpy), redondeado a n cifras significativas.
    Un vector (arreglo de una dimensión) se muestra como columna."""
    M = np.asarray(M, dtype=float)
    M = M.reshape(-1, 1) if M.ndim == 1 else M
    S = sp.Matrix(M.tolist())
    tex_nombre = tex if tex is not None else _a_tex(nombre)
    _emite(_Formula(f"{tex_nombre} \\approx {_tex_decimales(S, n)}\\qquad ({S.rows}\\times {S.cols})",
                    f"{nombre}  ({dims(S)}) ≈\n{_plano(S, lambda e: cifras(e, n))}"))
