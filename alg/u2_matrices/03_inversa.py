# ID: ALG-U2-NB03
# Notebook: alg/u2_matrices.ipynb · sección 2.3 propiedades del determinante y matriz inversa
# Repositorio: alg/u2_matrices/03_inversa.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Propiedades y matriz inversa
#
# **Propiedades del determinante** ($A$ y $B$ de $n\times n$, $k$ un número):
#
# - $\det(AB)=\det A\cdot\det B$ y $\det(A^T)=\det A$.
# - $\det(kA)=k^n\det A$: el factor $k$ sale una vez **por renglón**, no una sola vez.
# - Operaciones de renglón: intercambiar dos renglones cambia el signo; multiplicar un renglón por $k$ multiplica el determinante por $k$; sumar a un renglón un múltiplo de otro no lo cambia.
#
# **Inversa.** $A^{-1}$ es la matriz que cumple $AA^{-1}=A^{-1}A=I$. Existe **si y solo si** $\det A\neq 0$. Para una $2\times 2$:
#
# $$\begin{bmatrix}a&b\\c&d\end{bmatrix}^{-1}=\frac{1}{ad-bc}\begin{bmatrix}d&-b\\-c&a\end{bmatrix}.$$
#
# En general, con la matriz de cofactores $C=[C_{ij}]$ y la **adjunta** $\operatorname{adj}(A)=C^T$ (la transpuesta de la matriz de cofactores):
#
# $$A^{-1}=\frac{\operatorname{adj}(A)}{\det A},\qquad (AB)^{-1}=B^{-1}A^{-1}.$$
#
# Igual que con la transpuesta, la inversa de un producto **invierte el orden**. Si $y=Ax$, entonces $x=A^{-1}y$. Después de calcular una inversa a mano, comprueba siempre que $AA^{-1}=I$: es la forma más rápida de atrapar un signo o una transposición olvidada.

# %%
def matriz_cofactores(A):
    """C = [C_ij], con C_ij = (-1)^(i+j) M_ij."""
    A = como_matriz(A, "A")
    _cuadrada(A, "la matriz de cofactores solo existe para matrices cuadradas")
    n = A.rows
    if n == 1:
        return sp.Matrix([[1]])                     # convención: adj de una 1×1 es [1]
    return sp.Matrix(n, n, lambda i, j: cofactor(A, int(i) + 1, int(j) + 1))


def adjunta(A):
    """adj(A) = Cᵀ, la transpuesta de la matriz de cofactores."""
    return matriz_cofactores(A).T


def inversa(A):
    """(A⁻¹, det A) exactas por A⁻¹ = adj(A)/det A (hasta 7×7; más grande, con A.inv() de sympy); ValueError si det A = 0."""
    A = como_matriz(A, "A")
    _cuadrada(A, "la inversa solo existe para matrices cuadradas")
    d = A.det()
    if d == 0:
        raise ValueError("det A = 0: la matriz no tiene inversa")
    if A.rows > _MAX_COFACTORES:                      # n² cofactores de (n-1)×(n-1) ya no valen la pena
        return A.inv(), d
    return adjunta(A) / d, d


def es_identidad(M):
    """True si M es la matriz identidad."""
    M = como_matriz(M)
    return M.rows == M.cols and M == sp.eye(M.rows)


def vector_nulo(A):
    """Un x ≠ 0 con A x = 0 (existe cuando det A = 0), con entradas enteras."""
    A = como_matriz(A, "A")
    base = A.nullspace()
    if not base:
        raise ValueError("A x = 0 solo tiene la solución x = 0 (det A ≠ 0)")
    v = base[0]
    v = v * math.lcm(*[int(e.q) for e in v])         # sin fracciones
    v = v / math.gcd(*[int(e.p) for e in v])         # enteros lo más chicos posible
    primero = next(e for e in v if e != 0)
    return v if primero > 0 else -v                  # la primera entrada distinta de cero, positiva


def despeja(A, y):
    """x = A⁻¹ y, la solución exacta de y = A x."""
    A, y = como_matriz(A, "A"), vector(y, "y")
    _cuadrada(A, "para despejar x con A⁻¹, A debe ser cuadrada")
    if y.rows != A.rows:
        raise ValueError(f"y tiene {_entradas(y.rows)} y A tiene {A.rows} renglones: deben coincidir")
    try:
        Ainv, _ = inversa(A)
    except ValueError as err:
        raise ValueError(f"{err}, así que y = A x no tiene solución única "
                         "(puede no tener ninguna o tener infinitas; lo verás en la unidad de sistemas)") from None
    return Ainv * y


def calculadora_inversa(A_txt, n_cifras):
    try:
        A = matriz(A_txt, "A")
        _cuadrada(A, "la inversa solo existe para matrices cuadradas")
    except ValueError as err:
        print("Revisa la entrada:", err); return
    muestra("A", A, n_cifras)
    d = A.det()
    muestra("det A", d, n_cifras)
    if d == 0:
        print("det A = 0: A no tiene inversa. Hay un vector x ≠ 0 que A manda al vector cero, igual que manda 0 a 0;")
        print("como dos entradas distintas dan la misma salida, no hay forma de deshacer A. Por ejemplo:")
        x = vector_nulo(A)
        muestra("x", x)
        muestra("A x", A * x)
        return
    if A.rows > _MAX_COFACTORES:
        print(f"Con una {dims(A)}, la matriz de cofactores pide {A.rows ** 2} determinantes de "
              f"{A.rows - 1}×{A.rows - 1}; aquí se muestran hasta {_MAX_COFACTORES}×{_MAX_COFACTORES}.")
        print("La inversa se calcula con A.inv() de sympy (eliminación, también exacta).")
    else:
        print("Matriz de cofactores, C_ij = (-1)^(i+j)·M_ij:")
        muestra("C", matriz_cofactores(A))
        print("Adjunta: la transpuesta de C.")
        muestra("adj(A)", adjunta(A))
    Ainv, _ = inversa(A)
    if A.rows <= _MAX_COFACTORES:
        print("Inversa: la adjunta dividida entre det A.")
    muestra("A⁻¹", Ainv, n_cifras)
    P = A * Ainv
    muestra("A·A⁻¹", P)
    print("Verificación: A·A⁻¹ = I" if es_identidad(P) else "Verificación FALLÓ: A·A⁻¹ ≠ I")


widgets.interact(calculadora_inversa,
    A_txt=widgets.Textarea(value="2, 1, 0\n0, 2, 1\n1, 0, 2", description="A =", continuous_update=False,
                           layout=widgets.Layout(width="320px", height="110px")),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False));

# %% [markdown]
# **Despeja $x$ de $y=Ax$.** Si conoces la salida $y$ de un sistema lineal y la matriz $A$ que lo describe, la entrada que la produjo es $x=A^{-1}y$. La calculadora la obtiene exacta y comprueba que $Ax=y$.

# %%
def calculadora_despeje(A_txt, y_txt, n_cifras):
    try:
        A, y = matriz(A_txt, "A"), vector(y_txt, "y")
        x = despeja(A, y)
    except ValueError as err:     # det A = 0 no es un error de captura: el sistema no tiene solución única
        print("No hay solución única:" if str(err).startswith("det A = 0") else "Revisa la entrada:", err); return
    muestra("A", A, n_cifras)
    muestra("y", y, n_cifras)
    muestra("x = A⁻¹y", x, n_cifras)
    muestra("A x", A * x, n_cifras)
    print("Verificación: A x = y" if A * x == y else "Verificación FALLÓ: A x ≠ y")


widgets.interact(calculadora_despeje,
    A_txt=widgets.Textarea(value="2, 1, 0\n0, 2, 1\n1, 0, 2", description="A =", continuous_update=False,
                           layout=widgets.Layout(width="320px", height="110px")),
    y_txt=widgets.Text(value="11, 12, 7", description="y =", continuous_update=False),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False));

# %%
# Casos de prueba de la sección 3 (resultado conocido)
_S = matriz("2, 1, 0; 0, 2, 1; 1, 0, 2")
_A3, _B3 = matriz("1, 2; 3, 5"), matriz("2, 1; 1, 1")
PRUEBAS_3 = [
    ("inversa de [[2,1],[5,3]] = [[3,-1],[-5,2]]", lambda: inversa(matriz("2, 1; 5, 3"))[0] == sp.Matrix([[3, -1], [-5, 2]])),
    ("S = [[2,1,0],[0,2,1],[1,0,2]]: det 9 y adj S = [[4,-2,1],[1,4,-2],[-2,1,4]]",
     lambda: inversa(_S)[1] == 9 and adjunta(_S) == sp.Matrix([[4, -2, 1], [1, 4, -2], [-2, 1, 4]])),
    ("S⁻¹·(11, 12, 7) = (3, 5, 2), S·S⁻¹ = I y S no es la identidad",
     lambda: despeja(_S, "11, 12, 7") == vector("3, 5, 2") and es_identidad(_S * inversa(_S)[0]) and not es_identidad(_S)),
    ("[[1,2],[2,4]] (det 0) se rechaza", lambda: _mensaje(lambda: inversa(matriz("1, 2; 2, 4"))) ==
                                                "det A = 0: la matriz no tiene inversa"),
    ("con det 0, A x = 0 tiene un x ≠ 0, proporcional a (2, -1)",
     lambda: (lambda v: v[0] != 0 and 2 * v == v[0] * vector("2, -1") and matriz("1, 2; 2, 4") * v == sp.zeros(2, 1))(
             vector_nulo(matriz("1, 2; 2, 4")))),
    ("(AB)⁻¹ = B⁻¹A⁻¹ (y ≠ A⁻¹B⁻¹)",
     lambda: inversa(producto(_A3, _B3))[0] == producto(inversa(_B3)[0], inversa(_A3)[0])
             != producto(inversa(_A3)[0], inversa(_B3)[0])),
    ("det(2S) = 2³·det S = 72", lambda: det_cofactores(2 * _S)[0] == 8 * det_cofactores(_S)[0] == 72),
    ("det(AB) = det A·det B y det(Aᵀ) = det A",
     lambda: det_cofactores(producto(_A3, _B3))[0] == det_cofactores(_A3)[0] * det_cofactores(_B3)[0]
             and det_cofactores(transpuesta(_S))[0] == det_cofactores(_S)[0]),
    ("1×1: [5]⁻¹ = [1/5]; una no cuadrada se rechaza",
     lambda: inversa(matriz("5"))[0] == sp.Matrix([[sp.Rational(1, 5)]]) and _rechaza(lambda: inversa(matriz("1, 2, 3")))),
]
for nombre, prueba in PRUEBAS_3:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 3).** Calcula a mano la inversa de $\begin{bmatrix}4&3\\1&1\end{bmatrix}$ y escríbela como texto, renglón por renglón.

# %%
mi_inversa = None     # escribe tu resultado como texto, por ejemplo: "1, 0; 0, 1" (se aceptan fracciones como 2/3)

if mi_inversa is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    try:
        mia = matriz(str(mi_inversa), "tu inversa")
    except ValueError as err:
        print("Revisa tu matriz:", err)
    else:
        A = matriz("4, 3; 1, 1")
        ref, _ = inversa(A)
        print("La calculadora da:")
        muestra("A⁻¹", ref)
        print(compara(mia, ref, "Intercambia a y d, cambia el signo de b y c y divide entre det A = ad - bc."))
        if mia.shape == A.shape:
            print("Tu matriz por A da la identidad." if es_identidad(A * mia) else "Tu matriz por A no da la identidad.")
