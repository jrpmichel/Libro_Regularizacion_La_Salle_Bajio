# ID: DIF-U6-NB09
# Notebook: dif/u6_aplicaciones.ipynb · sección 6.9 civil
# Repositorio: dif/u6_aplicaciones/09_civil.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 9. Civil: $dV/dx=-w(x)$, $dM/dx=V(x)$; momento flexionante máximo
#
# Convención del libro: carga $w$ positiva hacia abajo y momento positivo con la concavidad hacia arriba. El momento tiene un extremo donde $V=0$; en apoyos simples, $M=0$. Escribe la posición como `x` (en metros).
#
# La calculadora recibe $M(x)$ (kN·m) y el claro $L$; obtiene $V$ y $w$, revisa los apoyos y encuentra el momento de mayor valor absoluto comparando los ceros de $V$ con los extremos.

# %%
def viga(M_txt, L):
    """Diccionario con V, w, M en los apoyos y (x, M) del momento de mayor valor absoluto en [0, L]."""
    M = parsear(M_txt)
    L = sp.nsimplify(L)
    if L <= 0:
        raise ValueError("el claro L debe ser positivo.")
    V = sp.diff(M, x)
    ceros = [c for c in sp.solveset(V, x, domain=sp.Interval(0, L))]
    cand = {c: M.subs(x, c) for c in [sp.Integer(0), L, *ceros]}
    xm = max(cand, key=lambda c: abs(float(cand[c])))
    return {"V": V, "w": sp.simplify(-sp.diff(V, x)), "apoyos": (M.subs(x, 0), M.subs(x, L)),
            "ceros_V": ceros, "x_max": xm, "M_max": sp.simplify(cand[xm])}


def calculadora_viga(M_txt, L, n_cifras):
    try:
        r = viga(M_txt, L)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print("V(x) =", r["V"], "kN    w(x) = -V' =", r["w"], "kN/m")
    print("M en los apoyos:", r["apoyos"], "(en apoyos simples debe ser 0)")
    print(f"V = 0 en x = {[cifras(c, n_cifras) for c in r['ceros_V']]}")
    print(f"momento máximo: {cifras(r['M_max'], n_cifras)} kN·m en x = {cifras(r['x_max'], n_cifras)} m")


widgets.interact(calculadora_viga,
    M_txt=widgets.Text(value="12*x - x**3/3", description="M(x) ="),
    L=widgets.FloatText(value=6.0, description="L (m)"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 9 (resultado conocido)
PRUEBAS_9 = [
    ("uniforme w = 10, L = 8: máx 80 en 4", lambda: (viga("40*x - 5*x**2", 8)["x_max"], viga("40*x - 5*x**2", 8)["M_max"]) == (4, 80)),
    ("uniforme: w = 10", lambda: viga("40*x - 5*x**2", 8)["w"] == 10),
    ("triangular: V = 12 - x², w = 2x", lambda: iguales(viga("12*x - x**3/3", 6)["V"], 12 - x**2) and iguales(viga("12*x - x**3/3", 6)["w"], 2*x)),
    ("triangular: máx 16√3 en √12", lambda: viga("12*x - x**3/3", 6)["x_max"] == sp.sqrt(12)
                                            and iguales(viga("12*x - x**3/3", 6)["M_max"], 16*sp.sqrt(3))),
    ("apoyos de la triangular: M = 0 en 0 y 6", lambda: viga("12*x - x**3/3", 6)["apoyos"] == (0, 0)),
    ("claro no positivo se rechaza", lambda: _rechaza(lambda: viga("x", 0))),
]
for nombre, prueba in PRUEBAS_9:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 9).** Carga uniforme con $w=6$ kN/m y $L=5$ m: calcula a mano $M_{\max}=wL^2/8$ y escribe el número. Para comprobar, escribe en la calculadora $M(x)=15x-3x^2$.

# %%
mi_M = None              # escribe un número

ref = viga("15*x - 3*x**2", 5)["M_max"]
if mi_M is None:                     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    print("La calculadora da:", ref, "=", cifras(ref, 4), "kN·m")
    print("coincide" if abs(mi_M - float(ref)) < 1e-6 else "NO coincide: el máximo está donde V = 15 - 6x = 0")
