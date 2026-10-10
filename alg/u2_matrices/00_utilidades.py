# ID: ALG-U2-NB00
# Notebook: alg/u2_matrices.ipynb · utilidades compartidas
# Repositorio: alg/u2_matrices/00_utilidades.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# # Álg. U2 · Matrices y determinantes
# **Libro de regularización en matemáticas para ingeniería · ULSA Bajío**
#
# Este notebook acompaña la unidad 2 del bloque de álgebra lineal. Tiene una sección por subtema, en el mismo orden que el libro:
#
# 1. Operaciones con matrices y transpuesta
# 2. Determinantes: Sarrus y cofactores
# 3. Propiedades y matriz inversa
# 4. Interpretación geométrica del determinante
# 5. Factorización LU
#
# **Cómo usarlo.** Ejecuta las celdas en orden (Entorno de ejecución → Ejecutar todas). Cada sección tiene una celda de teoría, una calculadora, sus casos de prueba y una celda **Contrasta**: antes de ejecutarla, resuelve a mano el caso que se indica y escribe tu resultado. La referencia de la calculadora aparece solo después de que escribes el tuyo.
#
# **Etiqueta del repositorio.** `ALG-U2` (bloque, unidad). Las celdas de código de cada sección vienen de `alg/u2_matrices/NN_*.py`; los fragmentos que el libro imprime, de `alg/u2_matrices/fragmentos/`.
#
# **Convenciones.** Punto decimal. Una matriz se escribe **por renglones**: los renglones se separan con `;` o con un salto de línea, y las entradas de cada renglón con comas o espacios. Por ejemplo, `1, 2; 3, 4` es la matriz $\begin{bmatrix}1&2\\3&4\end{bmatrix}$ (también se acepta la forma de Python `[[1, 2], [3, 4]]`). Puedes escribir fracciones como `2/3`; los decimales se convierten a fracciones exactas (`0.1` se lee como $\tfrac{1}{10}$), así que las cuentas internas no acumulan error de redondeo. Un vector se escribe como lista, `1, 2, 3`, y se trata como columna. Los índices $i,j$ empiezan en 1, como en el libro. Los decimales que se muestran se *redondean* (no se truncan) a las cifras significativas que elijas.

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
    return _desde_renglones(entradas, nombre)


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


def _cuadrada(A, mensaje="el determinante solo existe para matrices cuadradas"):
    if A.rows != A.cols:
        raise ValueError(f"{mensaje} (esta es {dims(A)})")


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
    e = q.adjusted()
    if e > 15 or e < -4:                      # muy grandes o muy chicos: notación científica, 1.234×10^-7
        return f"{format(q.scaleb(-e), 'f')}×10^{e}"
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


def _suma(valores):
    """'-16 - 5 + 0' en lugar de '-16 + -5 + 0'."""
    texto = str(valores[0])
    for v in valores[1:]:
        texto += f" - {-v}" if v < 0 else f" + {v}"
    return texto


def _en_linea(M):
    """Matriz en una línea con la misma notación de entrada: [2, 3; 5, -2]."""
    return "[" + "; ".join(", ".join(str(e) for e in M.row(i)) for i in range(M.rows)) + "]"


def _tupla(v, n=None):
    """Vector como (a, b, c); con n, además sus decimales redondeados a n cifras."""
    exacto = "(" + ", ".join(str(e) for e in v) + ")"
    if n is None or all(e.is_integer for e in v):
        return exacto
    return exacto + " ≈ (" + ", ".join(cifras(e, n) for e in v) + ")"


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
        ("−", "-"), ("det ", r"\det "), ("adj", r"\operatorname{adj}")]


def _a_tex(nombre):
    for a, b in _TEX:
        nombre = nombre.replace(a, b)
    return nombre


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
        if n is not None and not M.is_integer:
            tex += rf" \approx {cifras(M, n)}"
            texto += f" ≈ {cifras(M, n)}"
    display(_Formula(tex, texto))


def compara(mia, ref, pista=""):
    """'coinciden' o en qué entradas difiere tu matriz de la referencia, más una pista."""
    if mia.shape != ref.shape:
        return f"NO coinciden: tu matriz es {dims(mia)} y la de referencia es {dims(ref)}. {pista}".strip()
    malas = [f"({i + 1}, {j + 1})" for i in range(ref.rows) for j in range(ref.cols) if mia[i, j] != ref[i, j]]
    if not malas:
        return "coinciden"
    cuales = "la entrada" if len(malas) == 1 else "las entradas"
    return f"NO coinciden en {cuales} {', '.join(malas)} (renglón, columna). {pista}".strip()
