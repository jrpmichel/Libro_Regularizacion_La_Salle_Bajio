# ID: ALG-U5-NB00
# Notebook: alg/u5_ortogonalidad.ipynb · utilidades compartidas
# Repositorio: alg/u5_ortogonalidad/00_utilidades.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# # Álg. U5 · Ortogonalidad y mínimos cuadrados
# **Libro de regularización en matemáticas para ingeniería · ULSA Bajío**
#
# Este notebook acompaña la unidad 5 del bloque de álgebra lineal. Tiene una sección por subtema, en el mismo orden que el libro:
#
# 1. Producto interior, longitud y ortogonalidad
# 2. Proyección ortogonal sobre una recta y sobre un subespacio
# 3. Proceso de Gram-Schmidt
# 4. Mínimos cuadrados y ajuste de modelos lineales
#
# **Cómo usarlo.** Ejecuta las celdas en orden (Entorno de ejecución → Ejecutar todas). Cada sección tiene una celda de teoría, código que la recorre paso a paso, una calculadora, sus casos de prueba y una celda **Contrasta**: antes de ejecutarla, resuelve a mano el caso que se indica y escribe tu resultado. La referencia de la calculadora aparece solo después de que escribes el tuyo.
#
# **Etiqueta del repositorio.** `ALG-U5` (bloque, unidad). Las celdas de código de cada sección vienen de `alg/u5_ortogonalidad/NN_*.py`; los fragmentos que el libro imprime, de `alg/u5_ortogonalidad/fragmentos/`.
#
# **Convenciones.** En esta unidad hay dos maneras de escribir datos, y cada calculadora dice cuál usa.
#
# - **Una lista de vectores** (secciones 1, 2 y 3) se escribe **un vector por renglón**, o con los vectores separados por `;`; las componentes se separan con comas o espacios. **Cada renglón escrito es un vector**, y la calculadora lo guarda como **una columna** de la matriz $A=[\,\mathbf v_1\ \mathbf v_2\ \cdots\ \mathbf v_k\,]$, de $n\times k$ ($n$ componentes, $k$ vectores). Por ejemplo, `1, 0, 1; 1, 1, 0` son los vectores $\mathbf v_1=(1,0,1)$ y $\mathbf v_2=(1,1,0)$, y
#
#   $$A=\begin{bmatrix}1&1\\0&1\\1&0\end{bmatrix}.$$
#
#   Se aceptan hasta 6 vectores de hasta 6 componentes, todos con el mismo número de componentes. También se acepta `(1, 0, 1), (1, 1, 0)` o `{(1, 0, 1), (1, 1, 0)}`.
# - Un vector suelto ($\mathbf x$, la dirección $\mathbf u$ de una recta) se escribe como lista de números: `1, 2, 3`.
# - **Los datos de la sección 4** se escriben como listas de números separados por espacios o comas: `0 1 2 3 4` o `0, 1, 2, 3, 4`. La lista de $x$ y la de $y$ deben tener el mismo número de datos (hasta 500). Para ajustar un plano, cada renglón es un dato escrito como `x, y, z`.
# - Punto decimal. Puedes escribir fracciones como `2/3` (también `2 / 3`); los decimales se convierten a fracciones exactas (`0.1` se lee como $\tfrac{1}{10}$), así que los productos interiores, las proyecciones y las ecuaciones normales se calculan sin error de redondeo (en la sección 4, con hasta 12 datos y grado hasta 3; con más datos, más grado o logaritmos se calcula en decimales y la calculadora lo dice). Se aceptan los signos menos tipográficos (− y –) y un espacio después del signo (`- 1`).
# - Los vectores que escribes se nombran $\mathbf v_1, \mathbf v_2,\dots$ en el orden en que los escribes (en lo que imprimen las calculadoras, `v1`, `v2`, ...); los que da Gram-Schmidt, $\mathbf u_1,\mathbf u_2,\dots$ (ortogonales) y $\mathbf q_1,\mathbf q_2,\dots$ (de longitud 1). Vectores, renglones, columnas y datos se cuentan desde 1.
# - Un ángulo pide **siempre** su unidad, `"grados"` o `"radianes"`: no hay unidad por omisión.
# - En lo que muestran las calculadoras, el signo menos es siempre `-`, las raíces se escriben con `√`, la longitud con `‖ ‖` y el seno con `sen`. Los decimales que se muestran se *redondean* (no se truncan) a las cifras significativas que elijas; si una fracción exacta resulta demasiado larga, la calculadora muestra solo sus decimales y lo avisa.
# - Si usas las funciones desde código: a `lista_vectores` (y a las funciones que reciben vectores) le puedes pasar texto, una lista de vectores como `[(1, 0, 1), (1, 1, 0)]` o una `Matrix` de sympy; una `Matrix` se lee **por columnas** (cada columna es un vector), igual que la matriz $A$ que arma la calculadora. A las funciones de la sección 4 les puedes pasar texto o listas de números.

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
# Utilidades compartidas por todo el notebook (las mismas de la unidad 4): leer números y matrices exactos,
# redondear y mostrar
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
# Utilidades de esta unidad: listas de vectores (cada vector es una COLUMNA), eliminación exacta, núcleos,
# textos de combinaciones, gráficas, pruebas y lectura de datos
_MAX_VECTORES, _MAX_COMPONENTES = 6, 6        # hasta 6 vectores de hasta 6 componentes


def lista_vectores(texto, nombre="la lista de vectores"):
    """Lee una lista de vectores escrita UN VECTOR POR RENGLÓN o con los vectores separados por ';'
    ('1, 0, 1; 2, 1, 1') y devuelve la matriz A de n×k que tiene esos vectores como COLUMNAS.
    También acepta una lista de vectores ([(1, 0, 1), (2, 1, 1)] o ['1, 0, 1', '2, 1, 1']) y una Matrix de sympy,
    que se toma tal cual: sus columnas ya son los vectores."""
    if isinstance(texto, sp.MatrixBase):
        A = como_matriz(texto, nombre)
        lista = [list(A[:, j]) for j in range(A.cols)]
    elif isinstance(texto, str):
        lista, vacia = _trozos(texto)
        if vacia is not None:
            raise ValueError(f"en {nombre}, el vector {vacia} tiene una componente vacía: "
                             "revisa que no haya dos comas seguidas ni una coma al principio o al final")
    elif isinstance(texto, (list, tuple)):
        lista = []
        for i, v in enumerate(texto, 1):
            if isinstance(v, (str, list, tuple, np.ndarray, sp.MatrixBase)):
                lista.append(list(vector(v, f"el vector {i} de {nombre}")))
            else:
                raise ValueError(f"en {nombre}, escribe cada vector como una lista de números, por ejemplo: "
                                 "[(1, 0, 1), (2, 1, 1)]")
    else:
        raise ValueError(f"escribe {nombre} con un vector por renglón o separados por «;», por ejemplo: 1, 0, 1; 2, 1, 1")
    if not lista or not lista[0]:
        raise ValueError(f"escribe {nombre} con un vector por renglón o separados por «;», por ejemplo: 1, 0, 1; 2, 1, 1")
    k, n = len(lista), len(lista[0])
    for i, v in enumerate(lista, 1):
        if len(v) != n:
            raise ValueError(f"en {nombre}, el vector {i} tiene {_componentes(len(v))} y el vector 1 tiene {n}: "
                             f"todos los vectores deben tener el mismo número de componentes{pista_coma(texto)}")
    if k > _MAX_VECTORES:
        raise ValueError(f"{nombre} tiene {k} vectores; esta calculadora acepta hasta {_MAX_VECTORES} "
                         "(recuerda: cada renglón escrito es un vector)")
    if n > _MAX_COMPONENTES:
        raise ValueError(f"en {nombre}, los vectores tienen {n} componentes; esta calculadora acepta hasta "
                         f"{_MAX_COMPONENTES}{pista_coma(texto)}")
    columnas = []
    for i, v in enumerate(lista, 1):
        col = []
        for j, e in enumerate(v, 1):
            try:
                col.append(numero(e))
            except ValueError as err:
                raise ValueError(f"en {nombre}, vector {i}, componente {j}: {err}") from None
        columnas.append(col)
    return sp.Matrix(columnas).T


def _es_cero(e):
    """True si e vale 0. Un racional se compara directo; una expresión exacta (con raíces, senos o cosenos) se
    simplifica y, si sympy no la reduce, se evalúa con 50 cifras."""
    e = sp.sympify(e)
    if e.is_Rational:
        return e == 0
    v = sp.simplify(e)
    if v == 0:
        return True
    if v.is_number:
        return bool(abs(sp.N(v, 50)) < sp.Float("1e-40", 50))
    return False


def _simplifica(M):
    """Simplifica las entradas que no son racionales (las racionales se quedan como están)."""
    return M.applyfunc(lambda e: e if e.is_Rational else sp.simplify(e))


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
    """Eliminación de Gauss y después de Jordan, exacta, sobre una copia de M (adaptada de la unidad 3).

    Los pivotes se buscan en las primeras n_a columnas; se intercambian renglones solo si el pivote es 0 y abajo hay
    una entrada distinta de cero (como a mano). Con con_b, al final se busca un pivote en la columna n_a (la de b):
    si aparece, hay un renglón 0 = c. Funciona también con entradas exactas no racionales (senos y cosenos de una
    rotación): para decidir si una entrada es cero, primero se simplifica.
    Devuelve un diccionario: E (escalonada), R (reducida), pivotes (columnas con pivote, desde 0), libres, r (rango),
    r_aum (rango con la columna b), imposible ((renglón, c) o None) y pasos, una lista de (tipo, texto, matriz)."""
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
                motivo = "ya no quedan renglones abajo" if f >= m else f"de R{f + 1} hacia abajo todo es cero"
                pasos.append(("nota", f"la columna {c + 1} no tiene pivote: {motivo}", None))
            continue
        if k != f:
            pasos.append(("nota", f"el pivote de la columna {c + 1} es cero; se intercambian R{f + 1} y R{k + 1}", None))
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
            pasos.append(("nota", f"R{f + 1} quedó como 0 = {E[f, c]}: ninguna combinación lo cumple", None))
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
    if len(pasos) == hechas:
        pasos.append(("nota", "la forma escalonada ya era reducida: no hace falta ninguna operación", None))
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


def _primitivo(v):
    """El vector racional v (no nulo) escalado a enteros sin factor común y con su primera entrada distinta de cero
    positiva: (1/2, -1, 0) -> (1, -2, 0); (-2, -1, 1) -> (2, 1, -1)."""
    v = sp.Matrix(v)
    m = math.lcm(*[int(sp.Rational(e).q) for e in v])
    w = [int(e * m) for e in v]
    g = math.gcd(*w)
    signo = 1 if next(e for e in w if e != 0) > 0 else -1
    return sp.Matrix([signo * e // g for e in w])


def _termino(c, nombre, primero):
    """'2·v1', '- 1·v3', '+ (1/2)·v2': el coeficiente siempre se escribe, también el 1 y el 0. Si la fracción exacta
    es muy larga, se escribe redondeada: '+ (≈1.234×10^5)·v2'."""
    c = sp.sympify(c)
    mag = abs(c)
    coef = _txt(mag) if mag.is_integer else f"({_txt(mag)})"
    if len(coef) > 30:
        coef = f"(≈{cifras(mag, 4)})"
    negativo = bool(c.is_real and c < 0)
    if primero:
        return f"{'-' if negativo else ''}{coef}·{nombre}"
    return f"{'-' if negativo else '+'} {coef}·{nombre}"


def _subespacio(dim, n):
    """Nombre del subespacio según su dimensión: 'una recta por el origen', 'un plano por el origen', 'todo R³'..."""
    if dim == 0:
        return "{0} (solo el origen)"
    if dim == n:
        return f"todo {_R(n)}"
    if dim == 1:
        return "una recta por el origen"
    if dim == 2 and n == 3:
        return "un plano por el origen"
    return f"un subespacio de dimensión {dim} de {_R(n)} (pasa por el origen)"


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
    """Minúsculas, sin acentos ni signos de puntuación y con un solo espacio entre palabras: '¿Dependientes?' -> 'dependientes'."""
    import unicodedata
    t = unicodedata.normalize("NFD", str(texto).lower())
    t = "".join(c for c in t if unicodedata.category(c) != "Mn")
    return " ".join(re.sub(r"[^a-z0-9ñ ]+", " ", t).split())


# Gráficas: un color, un estilo de línea y un marcador por vector, en orden fijo, para que se lean también en blanco y negro
_COLORES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300"]
_ESTILOS = ["-", "--", ":", "-.", (0, (5, 1, 1, 1.5, 1, 1.5)), (0, (1, 4))]          # continuo, discontinuo, punteado...
_MARCAS = ["o", "s", "^", "D", "v", "P"]


def _flecha(ax, origen, punta, color, estilo, ancho, texto=None, desplazamiento=(6, 6)):
    """Flecha en el plano de origen a punta; con texto, su rótulo junto a la punta (sobre un fondo blanco)."""
    if np.hypot(*(np.asarray(punta, float) - np.asarray(origen, float))) > 1e-12:
        ax.annotate("", xy=punta, xytext=origen,
                    arrowprops=dict(arrowstyle="-|>", color=color, lw=ancho, linestyle=estilo, shrinkA=0, shrinkB=0))
    if texto:
        ax.annotate(texto, xy=punta, xytext=desplazamiento, textcoords="offset points", fontsize=11,
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


# --- utilidades propias de esta unidad ------------------------------------------------------------
def _punto(u, v):
    """Producto interior exacto de dos vectores columna: u·v = uᵀv."""
    return (u.T * v)[0]


def _es_nulo(v):
    """True si v es el vector cero."""
    return all(e == 0 for e in v)


def _valor(x, n=4):
    """Un valor exacto como texto, con sus decimales si no es entero: '2', '1/3 ≈ 0.3333', '√6 ≈ 2.449'. Si la forma
    exacta es muy larga, solo los decimales: '≈ 1.234 (aviso)'."""
    x = sp.sympify(x)
    if len(_txt(x)) > 28:
        return f"≈ {cifras(x, n)} ({_LARGA})"
    return _txt(x) + ("" if x.is_integer else f" ≈ {cifras(x, n)}")


def _eq(nombre, texto):
    """'nombre = texto', o 'nombre ≈ ...' si el texto ya es aproximado (empieza con '≈')."""
    return f"{nombre} {texto}" if texto.startswith("≈") else f"{nombre} = {texto}"


def _cuenta(k, singular, plural):
    """'1 dato', '5 datos': el número con el sustantivo en singular o en plural."""
    return f"{k} {singular if k == 1 else plural}"


def _cociente(num, den, n=4):
    """'(3/2)/(1/2) = 3': el cociente con paréntesis si el numerador o el denominador no es entero, y su valor.
    Si el cociente ya es el valor reducido ('3/2'), o si es muy largo, solo el valor."""
    num, den = sp.sympify(num), sp.sympify(den)
    parte = lambda e: _txt(e) if e.is_Integer else f"({_txt(e)})"
    frac = f"{parte(num)}/{parte(den)}"
    if frac == _txt(num / den) or len(frac) > 40:
        return _valor(num / den, n)
    return f"{frac} = {_valor(num / den, n)}"


def _productos_nulos(nombres):
    """'v1·v2 = 0', 'v1·v2 = v1·v3 = v2·v3 = 0'; con cuatro o más vectores, 'el producto de cada par de vectores
    distintos es 0'."""
    if len(nombres) > 3:
        return "el producto de cada par de vectores distintos es 0"
    return " = ".join(f"{a}·{b}" for i, a in enumerate(nombres) for b in nombres[i + 1:]) + " = 0"


def _aviso_coma(*textos):
    """Aviso (aunque no haya error) si la entrada tiene comas entre dígitos y números separados solo por espacios,
    que es como se escriben los decimales con coma: la calculadora leyó cada coma como separador."""
    if any(pista_coma(t) for t in textos):
        print("Aviso: hay comas entre dígitos (1,5) y números separados solo por espacios. Si escribiste decimales "
              "con coma, usa punto decimal (1.5, no 1,5): la calculadora leyó cada coma como separador.")


def _nombra_vectores(A, letra="v"):
    """'Vectores: v1 = (...), v2 = (...)  (2 vectores de R⁴)'."""
    print(("Vector: " if A.cols == 1 else "Vectores: ") +
          ", ".join(f"{letra}{j + 1} = {_tupla(A[:, j])}" for j in range(A.cols)) +
          f"  ({A.cols} {'vector' if A.cols == 1 else 'vectores'} de {_R(A.rows)})")


def _tex_A(k):
    """LaTeX de A = [v1 v2 ... vk] con los vectores en negritas y separados."""
    return r"A=[\," + r"\ ".join(rf"\mathbf v_{{{j}}}" for j in range(1, k + 1)) + r"\,]"


def _nombre_A(k):
    return f"A = [{' '.join(f'v{j}' for j in range(1, k + 1))}]"


_MAX_DATOS = 500                 # datos de la sección 4: hasta 500 por lista


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
    exactos: '0 1 2 3' -> [0, 1, 2, 3]; '10.0, 12.1' -> [10, 121/10]. Acepta también una lista o un arreglo de
    números. Hasta 500 datos; los mensajes de error dicen qué dato está mal. Con con_texto=True devuelve también
    los datos como se escribieron ('10.0', '12.1'), para mostrarlos igual."""
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
