# ID: ALG-U2-NB05
# Notebook: alg/u2_matrices.ipynb · sección 2.5 factorización LU
# Repositorio: alg/u2_matrices/05_lu.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 5. Factorización LU
#
# En el libro esta sección es una ampliación opcional. La eliminación gaussiana, sin intercambiar renglones, escribe una matriz cuadrada como $A=LU$ (factorización de Doolittle):
#
# - $U$ es **triangular superior**: es la matriz que queda al final de la eliminación.
# - $L$ es **triangular inferior con unos en la diagonal** y guarda los multiplicadores: en el paso $k$, al renglón $i$ le restas $\ell_{ik}$ veces el renglón del pivote, con $\ell_{ik}=\dfrac{\text{entrada }(i,k)\text{ en ese paso}}{\text{pivote }u_{kk}}$, y ese $\ell_{ik}$ queda en la posición $(i,k)$ de $L$.
# - Como $\det L=1$, el determinante sale gratis: $\det A=\det U=u_{11}u_{22}\cdots u_{nn}$, con unas $n^3/3$ multiplicaciones en lugar de las $\approx 1.7\,n!$ de los cofactores.
#
# Para resolver $A\mathbf x=\mathbf b$ se hacen dos sustituciones:
#
# $$L\mathbf y=\mathbf b\ \ (\text{hacia adelante: } y_1, y_2, \dots),\qquad U\mathbf x=\mathbf y\ \ (\text{hacia atrás: } x_n, x_{n-1}, \dots).$$
#
# La ventaja aparece cuando la misma $A$ se usa con muchos $\mathbf b$: se factoriza una sola vez y cada $\mathbf b$ nuevo cuesta solo las dos sustituciones.
#
# Si en algún paso el pivote $u_{kk}$ vale 0 y hay entradas distintas de cero debajo de él, hay que **intercambiar renglones**; con la matriz de permutación $P$ que registra los intercambios queda $PA=LU$. Eso es lo que hacen `numpy` y `scipy` (que además intercambian renglones para usar el pivote más grande y reducir el error de redondeo). Esta calculadora no intercambia renglones: te avisa en qué paso haría falta.
#
# En la caja de $\mathbf b$ puedes escribir uno o varios lados derechos separados con `;`, por ejemplo `10, 20, 30; 0, 0, 10`.

# %%
def lu_sin_pivoteo(A):
    """(L, U) exactas con A = LU (Doolittle: L con unos en la diagonal, sin intercambiar renglones)."""
    A = como_matriz(A, "A")
    _cuadrada(A, "aquí la factorización LU es solo para matrices cuadradas")
    n = A.rows
    L, U = sp.eye(n), A.copy()
    for k in range(n - 1):
        if U[k, k] == 0:
            if any(U[i, k] != 0 for i in range(k + 1, n)):
                raise ValueError(f"pivote cero en el paso {k + 1}: hace falta intercambiar renglones (PA = LU)")
            continue                                   # debajo del pivote ya hay ceros: no hay nada que eliminar
        for i in range(k + 1, n):
            m = U[i, k] / U[k, k]
            L[i, k] = m
            U[i, :] = U[i, :] - m * U[k, :]
    return L, U


def det_triangular(T):
    """Determinante de una matriz triangular: el producto de su diagonal."""
    T = como_matriz(T, "U")
    _cuadrada(T)
    if not (T.is_upper or T.is_lower):
        raise ValueError("det_triangular necesita una matriz triangular")
    return sp.prod([T[k, k] for k in range(T.rows)])


def _revisa_triangular(T, b, letra, inferior):
    T, b = como_matriz(T, letra), vector(b, "b" if letra == "L" else "y")
    _cuadrada(T, f"{letra} debe ser cuadrada")
    if b.rows != T.rows:
        raise ValueError(f"el vector tiene {_entradas(b.rows)} y {letra} tiene {T.rows} renglones: deben coincidir")
    if inferior and not T.is_lower:
        raise ValueError("L debe ser triangular inferior (ceros arriba de la diagonal)")
    if not inferior and not T.is_upper:
        raise ValueError("U debe ser triangular superior (ceros debajo de la diagonal)")
    return T, b


def adelante(L, b):
    """Sustitución hacia adelante: resuelve L y = b (L triangular inferior)."""
    L, b = _revisa_triangular(L, b, "L", inferior=True)
    n = L.rows
    y = sp.zeros(n, 1)
    for i in range(n):
        if L[i, i] == 0:
            raise ValueError(f"l{_sub(i + 1, i + 1)} = 0: no se puede despejar y{_sub(i + 1)}")
        y[i] = (b[i] - sum((L[i, j] * y[j] for j in range(i)), sp.Integer(0))) / L[i, i]
    return y


def atras(U, y):
    """Sustitución hacia atrás: resuelve U x = y (U triangular superior)."""
    U, y = _revisa_triangular(U, y, "U", inferior=False)
    n = U.rows
    x = sp.zeros(n, 1)
    for i in reversed(range(n)):
        if U[i, i] == 0:
            raise ValueError(f"u{_sub(i + 1, i + 1)} = 0 en la diagonal de U: det A = 0 y el sistema no tiene solución única")
        x[i] = (y[i] - sum((U[i, j] * x[j] for j in range(i + 1, n)), sp.Integer(0))) / U[i, i]
    return x


def resuelve_lu(L, U, b):
    """(y, x): primero L y = b hacia adelante, luego U x = y hacia atrás."""
    y = adelante(L, b)
    return y, atras(U, y)


def lados_derechos(texto, n):
    """Lista de vectores b de n entradas: uno por renglón ('10, 20, 30; 0, 0, 10') o uno solo escrito en columna."""
    try:
        B = matriz(texto, "b")
    except ValueError as err:
        if str(err).startswith("escribe "):
            raise ValueError("escribe b como lista de números, por ejemplo: 1, 2, 3 (varios b se separan con ;)") from None
        raise
    if B.cols == 1 and B.rows == n and n > 1:
        return [B]                                     # un solo b escrito en columna: 10; 20; 30
    if B.cols != n:
        raise ValueError(f"cada b debe tener {_entradas(n)}, tantas como renglones de A (escribiste {_entradas(B.cols)})")
    return [B.row(k).T for k in range(B.rows)]


def _numpy_resuelve(A, b, n_cifras):
    if A.det() == 0:          # en punto flotante, solve puede "resolver" una matriz singular y dar números sin sentido
        return "numpy.linalg.solve no aplica: det A = 0 (en punto flotante podría devolver números sin sentido sin avisar)"
    try:
        xs = np.linalg.solve(np.array(A.tolist(), dtype=float), np.array(b.tolist(), dtype=float).ravel())
    except np.linalg.LinAlgError:
        return "numpy.linalg.solve: la matriz es singular (det A = 0)"
    xs = np.where(np.abs(xs) < 1e-12 * max(1.0, float(np.abs(xs).max())), 0.0, xs)   # residuos de redondeo como 1e-16 -> 0
    return "numpy.linalg.solve (con intercambio de renglones, en punto flotante): x ≈ (" + ", ".join(cifras(v, n_cifras) for v in xs) + ")"


def calculadora_lu(A_txt, b_txt, n_cifras):
    try:
        A = matriz(A_txt, "A")
        _cuadrada(A, "aquí la factorización LU es solo para matrices cuadradas")
        bs = lados_derechos(b_txt, A.rows)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    muestra("A", A, n_cifras)
    try:
        L, U = lu_sin_pivoteo(A)
    except ValueError as err:
        print("Sin intercambiar renglones no se puede:", err)
        if A.det() == 0:
            print("Además det A = 0: ni intercambiando renglones hay solución única.")
            return
        print("numpy sí resuelve el sistema porque intercambia renglones (PA = LU):")
        for k, b in enumerate(bs, 1):
            print(f"  b{_sub(k)} = {_tupla(b)}:", _numpy_resuelve(A, b, n_cifras))
        return
    muestra("L", L, n_cifras)
    muestra("U", U, n_cifras)
    print("Verificación: LU = A" if L * U == A else "Verificación FALLÓ: LU ≠ A")
    d, n = det_triangular(U), U.rows
    nombres = "·".join(f"u{_sub(k, k)}" for k in range(1, n + 1)) if n <= 3 else f"u₁₁·…·u{_sub(n, n)}"
    cuenta = f" = {'·'.join(_paren(U[k, k]) for k in range(n))}" if n > 1 else ""
    print(f"det A = {nombres}{cuenta} = {d}" + ("" if d.is_integer else f" ≈ {cifras(d, n_cifras)}"))
    for k, b in enumerate(bs, 1):
        print(f"\n— Lado derecho b{_sub(k)} = {_tupla(b)} —")
        try:
            y, x = resuelve_lu(L, U, b)
        except ValueError as err:
            print("No se puede despejar:", err)
            print(_numpy_resuelve(A, b, n_cifras))
            continue
        print("L y = b (hacia adelante):  y =", _tupla(y, n_cifras))
        print("U x = y (hacia atrás):     x =", _tupla(x, n_cifras))
        print("Verificación: A x = b" if A * x == b else "Verificación FALLÓ: A x ≠ b")
        print(_numpy_resuelve(A, b, n_cifras))


widgets.interact(calculadora_lu,
    A_txt=widgets.Textarea(value="4, -2, 0\n-2, 4, -2\n0, -2, 2", description="A =", continuous_update=False,
                           layout=widgets.Layout(width="320px", height="110px")),
    b_txt=widgets.Textarea(value="10, 20, 30; 0, 0, 10; 5, 5, 5", description="b =", continuous_update=False,
                           layout=widgets.Layout(width="320px", height="60px")),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False));

# %%
# Casos de prueba de la sección 5 (resultado conocido)
_K = matriz("4, -2, 0; -2, 4, -2; 0, -2, 2")
_LK, _UK = lu_sin_pivoteo(_K)
_A4 = matriz("1, 2, 1, 0; 2, 5, 3, 1; 1, 3, 4, 2; 0, 1, 2, 5")
PRUEBAS_5 = [
    ("L = [[1,0,0],[-1/2,1,0],[0,-2/3,1]] y U = [[4,-2,0],[0,3,-2],[0,0,2/3]]",
     lambda: _LK == matriz("1, 0, 0; -1/2, 1, 0; 0, -2/3, 1") and _UK == matriz("4, -2, 0; 0, 3, -2; 0, 0, 2/3")),
    ("det A desde U: 4·3·(2/3) = 8", lambda: det_triangular(_UK) == 8 == _K.det()),
    ("b = (10, 20, 30) -> x = (30, 55, 70)", lambda: resuelve_lu(_LK, _UK, "10, 20, 30")[1] == vector("30, 55, 70")),
    ("b = (0, 0, 10) -> x = (5, 10, 15)", lambda: resuelve_lu(_LK, _UK, "0, 0, 10")[1] == vector("5, 10, 15")),
    ("b = (5, 5, 5) -> x = (15/2, 25/2, 15)", lambda: resuelve_lu(_LK, _UK, "5, 5, 5")[1] == vector("15/2, 25/2, 15")),
    ("[[0,1],[1,0]] se rechaza: pivote cero en el paso 1",
     lambda: _mensaje(lambda: lu_sin_pivoteo(matriz("0, 1; 1, 0"))).startswith("pivote cero en el paso 1")),
    ("pivote cero a la mitad (paso 2 de una 3×3) también se rechaza",
     lambda: "paso 2" in _mensaje(lambda: lu_sin_pivoteo(matriz("1, 2, 3; 2, 4, 5; 1, 1, 1")))),
    ("LU = A para una 4×4 (L con unos en la diagonal, U triangular superior)",
     lambda: (lambda L, U: L * U == _A4 and L.is_lower and U.is_upper and all(L[k, k] == 1 for k in range(4)))(*lu_sin_pivoteo(_A4))),
    ("varios b separados con ';' y un b escrito en columna",
     lambda: len(lados_derechos("10, 20, 30; 0, 0, 10", 3)) == 2 and lados_derechos("10; 20; 30", 3)[0] == vector("10, 20, 30")),
]
for nombre, prueba in PRUEBAS_5:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 5).** Factoriza a mano $A=\begin{bmatrix}3&1\\6&5\end{bmatrix}=LU$ (sin intercambiar renglones) y escribe $L$ y $U$ como texto.

# %%
mi_L = None     # escribe L como texto, por ejemplo: "1, 0; 5, 1"
mi_U = None     # escribe U como texto, por ejemplo: "3, 2; 0, 7"

if mi_L is None or mi_U is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    try:
        mia_L, mia_U = matriz(str(mi_L), "tu L"), matriz(str(mi_U), "tu U")
    except ValueError as err:
        print("Revisa tu matriz:", err)
    else:
        ref_L, ref_U = lu_sin_pivoteo(matriz("3, 1; 6, 5"))
        print("La calculadora da:")
        muestra("L", ref_L)
        muestra("U", ref_U)
        print("L:", compara(mia_L, ref_L, "El multiplicador es ℓ₂₁ = a₂₁/a₁₁ y va debajo del 1 de la diagonal."))
        print("U:", compara(mia_U, ref_U, "U es lo que queda al restar ℓ₂₁ veces el renglón 1 al renglón 2."))
