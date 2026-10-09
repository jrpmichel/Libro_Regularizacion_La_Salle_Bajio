# ID: DIF-U6-NB03
# Notebook: dif/u6_aplicaciones.ipynb · sección 6.3 rapideces relacionadas
# Repositorio: dif/u6_aplicaciones/03_rapideces.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Rapideces de variación relacionadas
#
# Si $x(t)$ y $y(t)$ cumplen $F(x,y)=0$ en todo instante, al derivar respecto a $t$:
# $$F_x\,\frac{dx}{dt}+F_y\,\frac{dy}{dt}=0\quad\Longrightarrow\quad \frac{dy}{dt}=-\frac{F_x}{F_y}\,\frac{dx}{dt}.$$
# Pasos: dibujar y nombrar; escribir la relación que vale **siempre**; derivar; sustituir los valores del instante **después** de derivar; interpretar signo y unidades.
#
# La calculadora recibe la ecuación (con `=`), un punto de la curva y $dx/dt$, y devuelve $dy/dt$.

# %%
def rapidez(ecuacion, px, py, dxdt):
    """dy/dt en el punto (px, py) de la curva, conocida dx/dt."""
    if ecuacion.count("=") != 1:
        raise ValueError("escribe una ecuación con un solo signo =, por ejemplo: x**2 + y**2 = 25")
    izq, der = ecuacion.split("=")
    F = parsear(izq, ("x", "y")) - parsear(der, ("x", "y"))
    punto = {x: sp.nsimplify(px), y: sp.nsimplify(py)}
    if sp.simplify(F.subs(punto)) != 0:
        raise ValueError(f"el punto ({px}, {py}) no cumple la ecuación: no corresponde a ningún instante.")
    Fx, Fy = sp.diff(F, x).subs(punto), sp.diff(F, y).subs(punto)
    if sp.simplify(Fy) == 0:
        raise ValueError("el coeficiente de dy/dt se anula en ese punto: la relación no determina dy/dt.")
    return sp.simplify(-Fx / Fy * sp.nsimplify(dxdt))


def calculadora_rapidez(ecuacion, px, py, dxdt, n_cifras):
    try:
        v = rapidez(ecuacion, px, py, dxdt)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    print(f"dy/dt = {v} ≈ {cifras(v, n_cifras)}  ({n_cifras} cifras significativas)")
    print("signo:", "y crece" if v > 0 else "y decrece" if v < 0 else "y no cambia en ese instante")


widgets.interact(calculadora_rapidez,
    ecuacion=widgets.Text(value="x**2 + y**2 = 25", description="ecuación"),
    px=widgets.FloatText(value=3.0, description="x"), py=widgets.FloatText(value=4.0, description="y"),
    dxdt=widgets.FloatText(value=0.5, description="dx/dt"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 3 (resultado conocido). En la tolva, x es el volumen y y la altura: x = π y³/12.
PRUEBAS_3 = [
    ("escalera en (3, 4) con dx/dt = 0.5: -0.375", lambda: rapidez("x**2 + y**2 = 25", 3, 4, 0.5) == sp.Rational(-3, 8)),
    ("escalera en (4, 3) con dx/dt = 0.3: -0.4", lambda: rapidez("x**2 + y**2 = 25", 4, 3, 0.3) == sp.Rational(-2, 5)),
    ("tolva h = 2 con dV/dt = -0.4: -0.4/π",
     lambda: iguales(rapidez("x = pi*y**3/12", sp.pi*8/12, 2, -0.4), -sp.Rational(2, 5)/sp.pi)),
    ("tolva h = 1: cuatro veces más rápido",
     lambda: iguales(rapidez("x = pi*y**3/12", sp.pi/12, 1, -0.4) / rapidez("x = pi*y**3/12", sp.pi*8/12, 2, -0.4), 4)),
    ("punto fuera de la curva se rechaza", lambda: _rechaza(lambda: rapidez("x**2 + y**2 = 25", 1, 1, 1))),
    ("coeficiente nulo se avisa", lambda: _rechaza(lambda: rapidez("x**2 + y**2 = 25", 5, 0, 1))),
]
for nombre, prueba in PRUEBAS_3:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 3).** Escalera de 5 m: calcula a mano $dy/dt$ cuando la base está a 4 m de la pared y se aleja a 0.4 m/s. Escribe el número.

# %%
mi_valor = None          # escribe un número

ref = rapidez("x**2 + y**2 = 25", 4, 3, 0.4)
if mi_valor is None:                     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    print("La calculadora da:", ref, "≈", cifras(ref, 4), "m/s")
    print("coincide" if abs(mi_valor - float(ref)) < 5e-4 else "NO coincide: 2x x' + 2y y' = 0 con x = 4, y = 3")
