# ID: ALG-U3-NB04
# Notebook: alg/u3_sistemas.ipynb · sección 3.4 clasificación de sistemas y sistemas homogéneos
# Repositorio: alg/u3_sistemas/04_clasificacion.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 4. Clasificación y sistemas homogéneos
#
# Con dos rangos, el de $A$ y el de la aumentada, clasificas cualquier sistema sin terminar de resolverlo.
#
# **Teorema de Rouché-Frobenius.** Sea $A$ de $m\times n$, $r=\operatorname{rango}(A)$ y $r^*=\operatorname{rango}[A\,|\,\mathbf b]$. Siempre $r\le r^*\le r+1$, y:
#
# | Rangos | Clasificación | Soluciones |
# |---|---|---|
# | $r<r^*$ | incompatible | ninguna: al escalonar aparece un renglón $0=c$, con $c\neq0$ |
# | $r=r^*=n$ | compatible determinado | una sola |
# | $r=r^*<n$ | compatible indeterminado | infinitas, que dependen de $n-r$ parámetros |
#
# **Variables básicas y libres.** En la forma escalonada de un sistema compatible, las incógnitas de las columnas con pivote son *básicas* y las demás son *libres*. La **solución general** da cada básica en función de las libres, que se renombran como parámetros.
#
# *Convención de este notebook:* la primera variable libre (de izquierda a derecha) se llama $s$, la segunda $t$, y siguen $u$, $v$, $w$ (con más de cinco, $s_1, s_2,\dots$). Cada parámetro **es** el valor de su variable libre: si $y$ es libre, $y=s$.
#
# **Particular más homogénea.** Si $\mathbf x_p$ resuelve $A\mathbf x=\mathbf b$, todas las soluciones son $\mathbf x=\mathbf x_p+\mathbf x_h$, con $A\mathbf x_h=\mathbf 0$, porque la diferencia de dos soluciones cumple $A(\mathbf x-\mathbf x_p)=\mathbf b-\mathbf b=\mathbf 0$. Por eso un sistema compatible tiene una solución o infinitas, nunca dos ni tres. La calculadora escribe la solución general así: $\mathbf x=\mathbf x_p+s\,\mathbf d_1+t\,\mathbf d_2+\cdots$, con $\mathbf x_p$ la solución que tiene los parámetros en 0.
#
# **Sistemas homogéneos** ($\mathbf b=\mathbf 0$, es decir, $A\mathbf x=\mathbf 0$). Siempre son compatibles, porque $\mathbf x=\mathbf 0$ (la solución *trivial*) los cumple. Tienen soluciones distintas de la trivial si y solo si $r<n$; en particular, siempre que haya menos ecuaciones que incógnitas y, si $A$ es cuadrada, cuando $\det A=0$. Todas las soluciones se obtienen dando valores a los parámetros de la solución general: con todos en cero sale la trivial, y con cualquier otro valor, una no trivial. Como la particular es $\mathbf 0$, la solución general es $s\,\mathbf d_1+t\,\mathbf d_2+\cdots$; la calculadora muestra además cada dirección $\mathbf d_i$ escalada para que su primera componente distinta de cero sea 1.
#
# **Ojo.** No clasifiques contando ecuaciones: cuatro ecuaciones con tres incógnitas pueden tener solución única, y tres con tres pueden no tener ninguna. Tampoco uses `numpy.linalg.lstsq`: devuelve un vector aunque el sistema sea incompatible, sin avisar. La calculadora decide con los rangos exactos.

# %%
def _escala_primera(d):
    """El vector d dividido entre su primera componente distinta de cero (para que esa componente sea 1)."""
    primera = next(e for e in d if e != 0)
    return d / primera


def clasificar(Ab, explicar=False, pasos=False, n_cifras=4):
    """Clasifica [A | b] con el teorema de Rouché-Frobenius, con rangos exactos (sin lstsq ni pinv).

    Devuelve un diccionario con m, n, r, r_aum, tipo ('incompatible', 'compatible determinado' o 'compatible
    indeterminado'), num_parametros, las formas escalonada y reducida, la solución (única, o general con los
    parámetros s, t, u, ... en el orden de las variables libres), la particular, las direcciones, el renglón (R, c)
    de 0 = c si es incompatible (tomado de la forma escalonada) y, si es homogéneo, las direcciones escaladas para que
    su primera componente distinta de cero sea 1. Con explicar=True imprime el informe; con pasos=True, también
    cada operación."""
    Ab = aumentada(Ab)
    m, n = Ab.rows, Ab.cols - 1
    A, b = Ab[:, :n], Ab[:, n]
    el = _elimina(Ab, n, con_b=True)
    r, r_aum = el["r"], el["r_aum"]
    if r < r_aum:
        tipo = "incompatible"
    elif r == n:
        tipo = "compatible determinado"
    else:
        tipo = "compatible indeterminado"
    res = {"m": m, "n": n, "r": r, "r_aum": r_aum, "tipo": tipo, "num_parametros": 0 if tipo == "incompatible" else n - r,
           "escalonada": el["E"], "reducida": el["R"], "pivotes": [c + 1 for c in el["pivotes"]],
           "libres": [c + 1 for c in el["libres"]], "homogeneo": b == sp.zeros(m, 1), "solucion": None,
           "particular": None, "direcciones": [], "parametros": [], "renglon": el["imposible"], "escaladas": [],
           "residuo": None, "pasos": el["pasos"]}
    if tipo != "incompatible":
        x, xp, dirs, ps = solucion_general(el["R"], el["pivotes"], n)
        res.update(solucion=x, particular=xp, direcciones=dirs, parametros=ps, residuo=b - A * xp)
        if res["homogeneo"]:
            res["escaladas"] = [_escala_primera(d) for d in dirs]
    if explicar:
        _informe_clasificacion(Ab, res, n_cifras, pasos)
    return res


def _informe_clasificacion(Ab, res, n_cifras, pasos=False):
    m, n, r, r_aum = res["m"], res["n"], res["r"], res["r_aum"]
    A = Ab[:, :n]
    nombres = nombres_incognitas(n)
    _imprime_sistema(Ab)
    muestra_aumentada("[A | b]", Ab, n)
    if pasos:
        _imprime_pasos(res["pasos"], corte=n)
    muestra_aumentada("forma escalonada", res["escalonada"], n)
    muestra_aumentada("forma escalonada reducida", res["reducida"], n)
    print(f"r = rango(A) = {r},  r* = rango[A | b] = {r_aum},  n = {n} (incógnitas)")
    if res["tipo"] == "incompatible":
        f, c = res["renglon"]
        print(f"r < r*: el sistema es INCOMPATIBLE (teorema de Rouché-Frobenius). El renglón R{f} de la forma escalonada "
              f"dice 0 = {c}, y ningún valor de las incógnitas lo cumple.")
        if c != 1:
            print("En la forma escalonada reducida ese renglón queda como 0 = 1.")
        return
    if res["homogeneo"]:
        print("Sistema homogéneo (b = 0): siempre es compatible, porque x = 0 (la solución trivial) lo cumple.")
    if res["tipo"] == "compatible determinado":
        x = res["solucion"]
        decimales = "" if all(e.is_integer for e in x) else " ≈ " + _aprox(x, n_cifras)
        print(f"r = r* = n = {n}: COMPATIBLE DETERMINADO, una sola solución.")
        print("Solución: " + ", ".join(f"{v} = {e}" for v, e in zip(nombres, x)) + decimales)
        if res["homogeneo"]:
            print(f"Solo tiene la solución trivial (r = n = {n}): no hay soluciones no triviales.")
        print("Comprobación: b - A x =", _tupla(res["residuo"]))
        return
    k = res["num_parametros"]
    ps = res["parametros"]
    libres = [nombres[c - 1] for c in res["libres"]]
    basicas = [nombres[c - 1] for c in res["pivotes"]]
    print(f"r = r* = {r} < n = {n}: COMPATIBLE INDETERMINADO, infinitas soluciones con n - r = {k} "
          f"parámetro{'s' if k > 1 else ''}.")
    print(f"Variables básicas (columnas con pivote): {_coma(basicas) if basicas else 'ninguna'}. "
          f"Libres: {_coma(libres)} → " + ", ".join(f"{v} = {p}" for v, p in zip(libres, ps)) + ".")
    print("Solución general:")
    for v, e in zip(nombres, res["solucion"]):
        print(f"  {v} = {_lineal(e, ps)}")
    vec = "(" + ", ".join(nombres) + ")"
    suma = " + ".join(f"{p}·{_tupla(d)}" for p, d in zip(ps, res["direcciones"]))
    cero = res["particular"] == sp.zeros(n, 1)
    print(f"Forma vectorial: {vec} = " + (suma if cero else f"{_tupla(res['particular'])} + {suma}"))
    ceros = all(A * d == sp.zeros(m, 1) for d in res["direcciones"])
    print(f"Comprobación: con {'todos los parámetros en 0' if k > 1 else f'{ps[0]} = 0'}, la particular {_tupla(res['particular'])} da "
          f"b - A x = {_tupla(res['residuo'])}; " +
          ("cada dirección d cumple A d = 0 (es solución del homogéneo)." if ceros else "una dirección NO cumple A d = 0."))
    if res["homogeneo"]:
        print(f"Hay soluciones no triviales (r = {r} < n = {n}): todas las soluciones se obtienen dando valores a "
              f"{'los parámetros' if k > 1 else ps[0]}; con {'todos en cero' if k > 1 else str(ps[0]) + ' = 0'} sale la trivial.")
        escaladas = ", ".join(_tupla(v) for v in res["escaladas"])
        print(f"Con {'cada dirección escalada' if k > 1 else 'la dirección escalada'} para que su primera componente "
              f"distinta de cero sea 1: {escaladas}. Las soluciones son " +
              (f"los múltiplos de {escaladas}." if k == 1 else "las sumas de múltiplos de esos vectores."))


_CLASIFICACIONES = {      # respuestas que se aceptan, ya normalizadas (minúsculas, sin acentos ni puntuación)
    "incompatible": ["incompatible", "es incompatible", "sistema incompatible", "sin solucion", "no tiene solucion",
                     "no hay solucion", "ninguna solucion", "ninguna"],
    "compatible determinado": ["determinado", "determinada", "compatible determinado", "compatible determinada",
                               "unica", "solucion unica", "una solucion", "una sola solucion", "tiene solucion unica"],
    "compatible indeterminado": ["indeterminado", "indeterminada", "compatible indeterminado", "compatible indeterminada",
                                 "infinitas", "infinitas soluciones", "tiene infinitas soluciones"],
}


def _tipo_de_texto(texto):
    """La clasificación que escribió el alumno, comparando la frase completa (así 'no es única' no se confunde con
    'única'). ValueError si no es una de las respuestas que se aceptan."""
    t = normaliza_texto(texto)
    for tipo, frases in _CLASIFICACIONES.items():
        if t in frases:
            return tipo
    raise ValueError(f"«{texto}» no es una de las respuestas que reconoce esta celda: escribe incompatible "
                     "(o sin solución), determinado (o solución única) o indeterminado (o infinitas soluciones)")

# %% [markdown]
# **Recorrido paso a paso.** Tres sistemas, uno de cada tipo: $x+y+z=6$, $2x+2y+3z=15$, $x+3y+2z=13$; luego $x+2y-z=3$, $2x+4y+z=9$; luego $x+y+z=1$, $2x+2y+2z=3$, $x-y=0$. Después, el sistema homogéneo del balance de la combustión del propano, $a\,\mathrm{C_3H_8}+b\,\mathrm{O_2}\to c\,\mathrm{CO_2}+d\,\mathrm{H_2O}$: carbono $3a-c=0$, hidrógeno $8a-2d=0$, oxígeno $2b-2c-d=0$ (aquí $a,b,c,d$ son $x_1,\dots,x_4$).

# %%
for _texto in ("1, 1, 1, 6; 2, 2, 3, 15; 1, 3, 2, 13", "1, 2, -1, 3; 2, 4, 1, 9", "1, 1, 1, 1; 2, 2, 2, 3; 1, -1, 0, 0"):
    clasificar(_texto, explicar=True)
    print("\n" + "=" * 90 + "\n")
clasificar("3, 0, -1, 0, 0; 8, 0, 0, -2, 0; 0, 2, -2, -1, 0", explicar=True);

# %% [markdown]
# **Calculadora.** Escribe la matriz aumentada (la última columna es $\mathbf b$; para un sistema homogéneo, ceros). Prueba sistemas con más ecuaciones que incógnitas y al revés.

# %%
def calculadora_clasificacion(Ab_txt, ver_pasos, n_cifras):
    with bloque():                 # el resultado aparece completo, de una sola vez
        try:
            Ab = aumentada(Ab_txt)
            res = clasificar(Ab)
        except ValueError as err:
            print("Revisa la entrada:", err); return
        _informe_clasificacion(Ab, res, n_cifras, pasos=ver_pasos)


widgets.interact(calculadora_clasificacion,
    Ab_txt=widgets.Textarea(value="1, 1, 1, 1, 4\n1, -1, 1, -1, 0", description="[A | b] =", continuous_update=False,
                            layout=widgets.Layout(width="340px", height="120px")),
    ver_pasos=widgets.Checkbox(value=False, description="mostrar cada operación"),
    n_cifras=widgets.IntSlider(value=4, min=2, max=8, description="cifras", continuous_update=False));

# %%
# Casos de prueba de la sección 4 (resultado conocido)
_s = sp.Symbol("s")
_CL = lambda texto: clasificar(texto)
_cumple = lambda texto, r: sp.expand(aumentada(texto)[:, :-1] * r["solucion"] - aumentada(texto)[:, -1]) == sp.zeros(r["m"], 1)
PRUEBAS_4 = [
    ("1, 1, 1, 6; 2, 2, 3, 15; 1, 3, 2, 13 → compatible determinado, (1, 2, 3)",
     lambda: (lambda r: r["tipo"] == "compatible determinado" and r["solucion"] == vector("1, 2, 3"))(
         _CL("1, 1, 1, 6; 2, 2, 3, 15; 1, 3, 2, 13"))),
    ("1, 2, -1, 3; 2, 4, 1, 9 → indeterminado, 1 parámetro: x = 4 - 2s, y = s, z = 1",
     lambda: (lambda r: r["tipo"] == "compatible indeterminado" and r["num_parametros"] == 1
                        and r["solucion"] == sp.Matrix([4 - 2 * _s, _s, 1]))(_CL("1, 2, -1, 3; 2, 4, 1, 9"))),
    ("1, 1, 1, 1; 2, 2, 2, 3; 1, -1, 0, 0 → incompatible, r = 2, r* = 3",
     lambda: (lambda r: r["tipo"] == "incompatible" and (r["r"], r["r_aum"]) == (2, 3) and r["solucion"] is None)(
         _CL("1, 1, 1, 1; 2, 2, 2, 3; 1, -1, 0, 0"))),
    ("homogéneo 3, 0, -1, 0, 0; 8, 0, 0, -2, 0; 0, 2, -2, -1, 0 → soluciones múltiplos de (1, 5, 3, 4)",
     lambda: (lambda r: r["homogeneo"] and r["tipo"] == "compatible indeterminado"
                        and r["escaladas"] == [vector("1, 5, 3, 4")])(_CL("3, 0, -1, 0, 0; 8, 0, 0, -2, 0; 0, 2, -2, -1, 0"))),
    ("homogéneo 1, 2, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0 → solo la solución trivial",
     lambda: (lambda r: r["homogeneo"] and r["tipo"] == "compatible determinado" and r["escaladas"] == []
                        and r["solucion"] == sp.zeros(3, 1))(_CL("1, 2, 0, 0; 0, 1, 0, 0; 0, 0, 1, 0"))),
    ("1, 1, 1, 1, 4; 1, -1, 1, -1, 0 → indeterminado, 2 parámetros (s, t), y la solución general cumple el sistema",
     lambda: (lambda r: r["tipo"] == "compatible indeterminado" and r["num_parametros"] == 2
                        and [str(p) for p in r["parametros"]] == ["s", "t"]
                        and _cumple("1, 1, 1, 1, 4; 1, -1, 1, -1, 0", r))(_CL("1, 1, 1, 1, 4; 1, -1, 1, -1, 0"))),
    ("1, 0, 0, 1; 0, 1, 0, 2; 0, 0, 1, 3; 1, 1, 1, 6 (4 ecuaciones, 3 incógnitas) → determinado, (1, 2, 3)",
     lambda: (lambda r: r["tipo"] == "compatible determinado" and r["solucion"] == vector("1, 2, 3"))(
         _CL("1, 0, 0, 1; 0, 1, 0, 2; 0, 0, 1, 3; 1, 1, 1, 6"))),
    ("el mismo con 7 en lugar de 6 en el último renglón → incompatible, r = 3, r* = 4",
     lambda: (lambda r: r["tipo"] == "incompatible" and (r["r"], r["r_aum"]) == (3, 4))(
         _CL("1, 0, 0, 1; 0, 1, 0, 2; 0, 0, 1, 3; 1, 1, 1, 7"))),
    ("0, 0, 0; 0, 0, 0 (2 ecuaciones, 2 incógnitas, todo cero) → indeterminado con 2 parámetros, r = 0",
     lambda: (lambda r: r["tipo"] == "compatible indeterminado" and r["r"] == 0 and r["num_parametros"] == 2)(
         _CL("0, 0, 0; 0, 0, 0"))),
    ("0, 0, 5 → incompatible: el renglón dice 0 = 5",
     lambda: (lambda r: r["tipo"] == "incompatible" and r["renglon"] == (1, 5))(_CL("0, 0, 5"))),
]
corre_pruebas(PRUEBAS_4, 4)

# %% [markdown]
# **Contrasta (sección 4).** Clasifica a mano $x-y=2$, $-2x+2y=-4$ con los rangos. Escribe si es incompatible, determinado o indeterminado y, si es indeterminado, cuántos parámetros tiene la solución general.

# %%
mi_clasificacion = None     # escribe "incompatible", "determinado" o "indeterminado" (también "sin solución", "solución única" o "infinitas soluciones")
mis_parametros = None       # si escribiste "indeterminado": cuántos parámetros, por ejemplo 2

if mi_clasificacion is None:     # la referencia se muestra después de tu cálculo
    print("falta tu clasificación: escríbela arriba y vuelve a ejecutar la celda")
else:
    try:
        mio = _tipo_de_texto(mi_clasificacion)
        mis_p = None if mis_parametros is None else numero(mis_parametros, "tu número de parámetros")
    except ValueError as err:
        print("Revisa tu respuesta:", err)
    else:
        if mio == "compatible indeterminado" and mis_p is None:     # respuesta incompleta: todavía sin referencia
            print("Si es indeterminado, escribe también cuántos parámetros tiene la solución general (mis_parametros) "
                  "y vuelve a ejecutar la celda.")
        elif mio == "compatible indeterminado" and not (mis_p.is_integer and mis_p >= 1):
            print(f"El número de parámetros es un entero de 1 en adelante (escribiste {mis_p}): revísalo y vuelve a "
                  "ejecutar la celda.")
        else:
            ref = clasificar("1, -1, 2; -2, 2, -4")
            texto_ref = ref["tipo"] + (f" con {ref['num_parametros']} parámetro{'s' if ref['num_parametros'] > 1 else ''}"
                                       if ref["num_parametros"] else "")
            print(f"La calculadora da: {texto_ref} (r = {ref['r']}, r* = {ref['r_aum']}, n = {ref['n']}).")
            if mio != ref["tipo"]:
                print("NO coinciden: escalona la aumentada y compara r = rango(A) con r* = rango[A | b] y con n.")
            elif ref["num_parametros"] and mis_p != ref["num_parametros"]:
                print("La clasificación coincide, pero el número de parámetros no: es n - r.")
            else:
                print("coinciden")
            print("\nEl informe completo:")
            clasificar("1, -1, 2; -2, 2, -4", explicar=True, pasos=True)
