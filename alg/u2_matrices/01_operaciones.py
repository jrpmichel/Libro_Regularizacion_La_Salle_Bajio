# ID: ALG-U2-NB01
# Notebook: alg/u2_matrices.ipynb · sección 2.1 operaciones con matrices y transpuesta
# Repositorio: alg/u2_matrices/01_operaciones.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 1. Operaciones con matrices y transpuesta
#
# Una matriz $A$ de tamaño $m\times n$ tiene $m$ renglones y $n$ columnas; su entrada $a_{ij}$ está en el renglón $i$ y la columna $j$.
#
# **Suma y resta** (solo entre matrices del mismo tamaño): $(A\pm B)_{ij}=a_{ij}\pm b_{ij}$. **Múltiplo escalar:** $(kA)_{ij}=k\,a_{ij}$.
#
# **Producto.** Si $A$ es $m\times n$ y $B$ es $n\times p$, la entrada $c_{ij}$ de $C=AB$ combina el renglón $i$ de $A$ con la columna $j$ de $B$:
#
# $$c_{ij}=\sum_{k=1}^{n}a_{ik}\,b_{kj}=a_{i1}b_{1j}+a_{i2}b_{2j}+\cdots+a_{in}b_{nj},\qquad (m\times n)(n\times p)=m\times p.$$
#
# Las columnas de $A$ deben ser tantas como los renglones de $B$ (los dos $n$ de en medio); los números de los extremos dan el tamaño del resultado.
#
# - **Identidad** $I_n$: unos en la diagonal y ceros fuera de ella. Si $A$ es $m\times n$, $I_mA=AI_n=A$, como el 1 entre los números.
# - **Transpuesta:** $(A^T)_{ij}=a_{ji}$; los renglones pasan a ser columnas. Cumple $(A^T)^T=A$, $(A+B)^T=A^T+B^T$ y $(AB)^T=B^TA^T$: al transponer un producto, **el orden se invierte**.
# - En general $AB\neq BA$. A veces los dos productos existen y dan matrices distintas; a veces uno de ellos ni siquiera está definido.
#
# **Ejemplo: materiales y costo.** Una planta fabrica tres productos con tres materiales. La entrada $m_{ij}$ de $M$ es la cantidad del material $i$ que lleva una unidad del producto $j$, y la entrada $p_{jk}$ de $P$ es el número de unidades del producto $j$ que se fabrican en la semana $k$. Entonces $MP$ ($3\times 3$ por $3\times 2$, resultado $3\times 2$) da el material $i$ que se consume en la semana $k$, y con el vector $\mathbf c$ de costos por unidad de material, el renglón $\mathbf c^T(MP)$ da el costo de material de cada semana:
#
# $$M=\begin{bmatrix}0.8&2.5&4\\12&20&32\\0.2&0.6&0.9\end{bmatrix},\quad P=\begin{bmatrix}40&25\\10&15\\6&8\end{bmatrix},\quad MP=\begin{bmatrix}81&89.5\\872&856\\19.4&21.2\end{bmatrix},\quad \mathbf c=\begin{bmatrix}180\\0.5\\120\end{bmatrix},\quad \mathbf c^T(MP)=\begin{bmatrix}17\,344&19\,082\end{bmatrix}.$$
#
# Supuesto del ejemplo: cada cantidad de $M$ va en la unidad propia de su material y $\mathbf c$ está en pesos por esa unidad, así que $\mathbf c^T(MP)$ queda en pesos. Para reproducirlo en la calculadora escribe A = `0.8, 2.5, 4; 12, 20, 32; 0.2, 0.6, 0.9`, B = `40, 25; 10, 15; 6, 8` y elige la operación AB.

# %%
def suma(A, B, nombres=("A", "B")):
    """A + B exacta; ValueError si las dos matrices no tienen el mismo tamaño."""
    a, b = nombres
    A, B = como_matriz(A, a), como_matriz(B, b)
    if A.shape != B.shape:
        raise ValueError(f"{a} es {dims(A)} y {b} es {dims(B)}: para sumar o restar, las dos deben tener el mismo tamaño")
    return A + B


def producto(A, B, nombres=("A", "B")):
    """AB exacto; ValueError que explica el tamaño si las columnas de A no son tantas como los renglones de B."""
    a, b = nombres
    A, B = como_matriz(A, a), como_matriz(B, b)
    if A.cols != B.rows:
        raise ValueError(f"{a} es {dims(A)} y {b} es {dims(B)}: para {a}{b}, las columnas de {a} ({A.cols}) "
                         f"deben ser tantas como los renglones de {b} ({B.rows})")
    return A * B


def transpuesta(A):
    """Aᵀ: el renglón i de A pasa a ser la columna i."""
    return como_matriz(A).T


def entrada_del_producto(A, B, i=1, j=1, nombres=("A", "B")):
    """Texto con la cuenta de la entrada (i, j) de AB: renglón i de A por columna j de B."""
    C = producto(A, B, nombres)
    A, B = como_matriz(A), como_matriz(B)
    terminos = " + ".join(f"{_paren(A[i - 1, k])}·{_paren(B[k, j - 1])}" for k in range(A.cols))
    return f"c{_sub(i, j)} = (renglón {i} de {nombres[0]})·(columna {j} de {nombres[1]}) = {terminos} = {C[i - 1, j - 1]}"


def ab_y_ba(A, B):
    """(AB o None, BA o None, veredicto en palabras)."""
    A, B = como_matriz(A, "A"), como_matriz(B, "B")
    try:
        AB = producto(A, B)
    except ValueError as err:
        AB, por_que_ab = None, str(err)
    try:
        BA = producto(B, A, ("B", "A"))
    except ValueError as err:
        BA, por_que_ba = None, str(err)
    if AB is None and BA is None:
        return AB, BA, f"ni AB ni BA están definidos: A es {dims(A)} y B es {dims(B)}"
    if BA is None:
        return AB, BA, f"AB sí está definido ({dims(AB)}), pero BA no: {por_que_ba}"
    if AB is None:
        return AB, BA, f"BA sí está definido ({dims(BA)}), pero AB no: {por_que_ab}"
    if AB.shape != BA.shape:
        return AB, BA, f"AB es {dims(AB)} y BA es {dims(BA)}: ni siquiera tienen el mismo tamaño, así que AB ≠ BA"
    if AB == BA:
        return AB, BA, "AB = BA: estas dos matrices conmutan (es la excepción, no la regla)"
    return AB, BA, "AB ≠ BA: el orden de los factores sí importa"


def _muestra_producto(A, B, nombres, n):
    C = producto(A, B, nombres)
    muestra("".join(nombres), C, n)
    print("Por ejemplo,", entrada_del_producto(A, B, 1, 1, nombres))


_OPERACIONES = ["A + B", "A − B", "AB", "BA", "AB y BA", "Aᵀ", "(AB)ᵀ y BᵀAᵀ"]


def calculadora_operaciones(A_txt, B_txt, operacion, n_cifras):
    try:
        A = matriz(A_txt, "A")
        B = None if operacion == "Aᵀ" else matriz(B_txt, "B")      # Aᵀ no necesita a B
        muestra("A", A, n_cifras)
        if B is not None:
            muestra("B", B, n_cifras)
        if operacion == "A + B":
            muestra("A + B", suma(A, B), n_cifras)
        elif operacion == "A − B":
            muestra("A − B", suma(A, -B), n_cifras)
        elif operacion == "AB":
            _muestra_producto(A, B, ("A", "B"), n_cifras)
        elif operacion == "BA":
            _muestra_producto(B, A, ("B", "A"), n_cifras)
        elif operacion == "AB y BA":
            AB, BA, veredicto = ab_y_ba(A, B)
            if AB is not None:
                muestra("AB", AB, n_cifras)
            if BA is not None:
                muestra("BA", BA, n_cifras)
            print(veredicto)
        elif operacion == "Aᵀ":
            muestra("Aᵀ", transpuesta(A), n_cifras)
            print(f"A es {dims(A)} y Aᵀ es {dims(transpuesta(A))}: el renglón 1 de A es la columna 1 de Aᵀ.")
        else:
            izquierda = transpuesta(producto(A, B))
            derecha = producto(transpuesta(B), transpuesta(A), ("Bᵀ", "Aᵀ"))
            muestra("(AB)ᵀ", izquierda, n_cifras)
            muestra("BᵀAᵀ", derecha, n_cifras)
            print("(AB)ᵀ = BᵀAᵀ: coinciden, como dice la regla." if izquierda == derecha else "NO coinciden: revisa la entrada.")
            try:
                AtBt = producto(transpuesta(A), transpuesta(B), ("Aᵀ", "Bᵀ"))
                if AtBt.shape != izquierda.shape:
                    print(f"Ojo: AᵀBᵀ (sin invertir el orden) es {dims(AtBt)}, ni siquiera del tamaño de (AB)ᵀ.")
                elif AtBt != izquierda:
                    print("Ojo: AᵀBᵀ (sin invertir el orden) da otra matriz:")
                    muestra("AᵀBᵀ", AtBt, n_cifras)
                else:
                    print("Aquí AᵀBᵀ también coincide, pero es casualidad: en general hay que invertir el orden.")
            except ValueError as err:
                print(f"Ojo: AᵀBᵀ (sin invertir el orden) ni siquiera está definido: {err}")
    except ValueError as err:
        print("Revisa la entrada:", err)


widgets.interact(calculadora_operaciones,
    A_txt=widgets.Text(value="1, 2; 3, 4", description="A =", continuous_update=False),
    B_txt=widgets.Text(value="0, 1; 1, 0", description="B =", continuous_update=False),
    operacion=widgets.Dropdown(options=_OPERACIONES, value="AB y BA", description="operación"),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False));

# %%
# Casos de prueba de la sección 1 (resultado conocido)
_A1, _B1 = matriz("1, 2; 3, 4"), matriz("0, 1; 1, 0")
_A23, _B32 = matriz("2, 0, 1; 1, -1, 3"), matriz("1, 2; 0, 1; 4, -2")
_M, _P = matriz("0.8, 2.5, 4; 12, 20, 32; 0.2, 0.6, 0.9"), matriz("40, 25; 10, 15; 6, 8")
PRUEBAS_1 = [
    ("[[1,2],[3,4]]·[[0,1],[1,0]] = [[2,1],[4,3]] y la cuenta de c₁₂ es 1·1 + 2·0 = 1",
     lambda: producto(_A1, _B1) == sp.Matrix([[2, 1], [4, 3]])
             and entrada_del_producto(_A1, _B1, 1, 2).endswith("= 1·1 + 2·0 = 1")),
    ("BA = [[3,4],[1,2]] y BA ≠ AB", lambda: producto(_B1, _A1) == sp.Matrix([[3, 4], [1, 2]])
                                            and ab_y_ba(_A1, _B1)[2].startswith("AB ≠ BA")),
    ("(2×3)(3×2) da una 2×2", lambda: dims(producto(_A23, _B32)) == "2×2"),
    ("(2×3)(2×3) se rechaza y el mensaje explica el tamaño",
     lambda: _mensaje(lambda: producto(_A23, _A23)) == "A es 2×3 y B es 2×3: para AB, las columnas de A (3) "
                                                     "deben ser tantas como los renglones de B (2)"),
    ("(AB)ᵀ = BᵀAᵀ", lambda: transpuesta(producto(_A23, _B32)) == producto(transpuesta(_B32), transpuesta(_A23))),
    ("materiales: MP = [[81, 89.5], [872, 856], [19.4, 21.2]]",
     lambda: producto(_M, _P) == matriz("81, 89.5; 872, 856; 19.4, 21.2")),
    ("costo por semana cᵀ(MP) = (17344, 19082)",
     lambda: producto(transpuesta(vector("180, 0.5, 120")), producto(_M, _P)) == sp.Matrix([[17344, 19082]])),
    ("A + B = [[1,3],[4,4]] y A + B de tamaños distintos se rechaza",
     lambda: suma(_A1, _B1) == sp.Matrix([[1, 3], [4, 4]]) and _rechaza(lambda: suma(_A1, _A23))),
    ("lectura exacta: 0.1 -> 1/10 y 2/3 queda como 2/3",
     lambda: matriz("0.1, 2/3")[0, 0] == sp.Rational(1, 10) and matriz("0.1 2/3")[0, 1] == sp.Rational(2, 3)),
    ("renglones disparejos, texto vacío y letras se rechazan",
     lambda: "el renglón 2 tiene 1 entrada y el 1 tiene 2" in _mensaje(lambda: matriz("1, 2; 3"))
             and _rechaza(lambda: matriz("  ")) and _rechaza(lambda: matriz("1, x; 3, 4"))),
]
for nombre, prueba in PRUEBAS_1:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 1).** Calcula a mano $\begin{bmatrix}1&-1\\2&0\end{bmatrix}\begin{bmatrix}3&1\\1&2\end{bmatrix}$ y escribe el resultado como texto, renglón por renglón.

# %%
mi_AB = None     # escribe tu resultado como texto, por ejemplo: "1, 0; 0, 1"

if mi_AB is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    try:
        mia = matriz(str(mi_AB), "tu resultado")
    except ValueError as err:
        print("Revisa tu matriz:", err)
    else:
        ref = producto(matriz("1, -1; 2, 0"), matriz("3, 1; 1, 2"))
        print("La calculadora da:")
        muestra("AB", ref)
        print(compara(mia, ref, "Cada c_ij es el renglón i de la primera matriz por la columna j de la segunda; "
                                "por ejemplo, c₁₂ = a₁₁b₁₂ + a₁₂b₂₂."))
