# ID: INT-U6-NB07
# Notebook: int/u6_aplicaciones.ipynb · sección 6.7 trabajo, energía y eficiencia
# Repositorio: int/u6_aplicaciones/07_trabajo.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 6.7 Trabajo, energía y eficiencia
#
# $$W=\int_{x_0}^{x_1}F(x)\,dx\quad(\text{J}=\text{N·m}),\qquad \eta=\frac{W_{\text{útil}}}{E_{\text{entrada}}},\quad 0<\eta\le1,\qquad 1\ \text{kWh}=3.6\times10^{6}\ \text{J}.$$
#
# Para un resorte, $F=kx$ y $W=\tfrac12k\,(x_1^2-x_0^2)$. Al bombear un líquido, cada rebanada sube una distancia distinta: en un tanque cilíndrico de radio $r$ lleno hasta $H$, la rebanada a la altura $y$ pesa $\rho g\,\pi r^2\,dy$ y sube $H-y$, así que el integrando es $\rho g\,\pi r^2(H-y)$. Si una cuenta da $\eta>1$, hay un error en los datos o en las unidades, porque ninguna máquina entrega más energía de la que recibe.
#
# El ejemplo exacto calcula el trabajo para estirar un resorte de $x_0$ a $x_1$.

# %%
k_s, x0_s, x1_s = sp.symbols("k x_0 x_1", positive=True)
W_resorte = sp.integrate(k_s * x, (x, x0_s, x1_s))
print("W =", W_resorte, "   con k = 800 N/m, de 0.05 a 0.15 m:",
      W_resorte.subs({k_s: 800, x0_s: sp.Rational(1, 20), x1_s: sp.Rational(3, 20)}), "J")

# %% [markdown]
# **Calculadora de trabajo.** Escribe $F(x)$ en N (también sirve la variable $y$), los extremos en m y, si la conoces, la energía de entrada en J o kWh (0 si no hay dato). La calculadora integra la fuerza, da el trabajo en J y, con la energía de entrada, la eficiencia. El ejemplo es un malacate que sube 30 m una carga de 400 kg con un cable de 1.5 kg/m; a la altura $y$ quedan colgando $30-y$ metros de cable.

# %%
def trabajo(F, x0, x1, con_metodo=False):
    """W = ∫_{x0}^{x1} F dx en J (F en N, x en m), exacto con sympy si se puede. Con con_metodo=True devuelve
    (W, método)."""
    F = sp.sympify(F)
    x0, x1 = exacto(x0, "x0"), exacto(x1, "x1")
    if x0 == x1:
        raise ValueError("x0 y x1 son iguales: no hay desplazamiento, el trabajo es 0.")
    revisa_intervalo(F, x0, x1, "F")
    W, como = integra(F, x0, x1)
    return (W, como) if con_metodo else W


def a_joules(valor, unidad="J"):
    """Energía en J a partir de J o kWh."""
    factores = {"J": 1, "kWh": 3600000}
    if unidad not in factores:
        raise ValueError('la unidad es "J" o "kWh".')
    return valor * factores[unidad]


def eficiencia(W_util, W_entrada):
    """η = W_util / W_entrada. ValueError si no cumple 0 < η <= 1."""
    W_util, W_entrada = float(W_util), float(W_entrada)
    if W_entrada <= 0:
        raise ValueError("la energía de entrada debe ser positiva.")
    if W_util <= 0:
        raise ValueError("el trabajo útil debe ser positivo para hablar de eficiencia.")
    eta = W_util / W_entrada
    if eta > 1:
        raise ValueError(f"saldría η = {cifras(eta, 4)} > 1: la máquina entregaría más energía de la que recibe, lo que "
                         "viola la conservación de la energía. Revisa los datos y las unidades (1 kWh = 3.6 MJ).")
    return eta


def calculadora_trabajo(F_txt, x0, x1, E_in, unidad):
    try:
        F = parsear(F_txt)
        W, como = trabajo(F, x0, x1, con_metodo=True)
    except ERRORES as err:
        print("Revisa la entrada:", explica(err)); return
    print(f"({CIFRAS} cifras significativas, redondeadas; {COMO[como]})")
    print(f"W = {muestra(W)} J" + (f" = {cifras(float(W) / 1000, CIFRAS)} kJ" if abs(float(W)) >= 1000 else "")
          + ("   (negativo porque x1 < x0: recorres el intervalo hacia atrás)" if x1 < x0 else ""))
    if E_in:
        try:
            E = a_joules(E_in, unidad)
            eta = eficiencia(W, E)
        except ValueError as err:
            print("Con la energía de entrada:", err)
        else:
            print(f"energía de entrada = {cifras(E, CIFRAS)} J   →   η = {cifras(eta, 4)} = {cifras(100 * eta, 3)} %"
                  "  (η con 4 cifras, porcentaje con 3)")
    try:
        xs = np.linspace(float(x0), float(x1), 400)
        ys = para_graficar(F, xs)
        fig, ax = plt.subplots(figsize=(5.5, 3))
        ax.fill_between(xs, 0, ys, color="0.85"); ax.plot(xs, ys, color="black")
        rotula(ax, xs[len(xs) // 2], ys[len(ys) // 2], "F(x)", xytext=(4, 8))
        ax.set_xlabel("x (m)"); ax.set_ylabel("F (N)"); ax.set_ylim(bottom=min(0, float(np.nanmin(ys))))
        ax.set_title("el área gris es el trabajo W", fontsize=9); plt.show()
    except ERRORES:
        print("(no pude dibujar la gráfica con estos datos)")


widgets.interact(calculadora_trabajo,
    F_txt=widgets.Text(value="400*9.81 + 1.5*9.81*(30 - y)", description="F =", continuous_update=False),
    x0=widgets.FloatText(value=0, description="x0 (m)"),
    x1=widgets.FloatText(value=30, description="x1 (m)"),
    E_in=widgets.FloatText(value=0.05, description="entrada"),
    unidad=widgets.Dropdown(options=["J", "kWh"], value="kWh", description="unidad"));

# %%
# Casos de prueba de la sección 6.7 (resultado conocido)
_g = sp.Rational("9.81")
_malacate = 400 * _g + sp.Rational("1.5") * _g * (30 - x)

PRUEBAS_7 = [
    ("resorte k = 200 N/m de 0 a 0.1 m: 1 J", lambda: trabajo(200 * x, 0, 0.1) == 1),
    ("resorte k = 800 N/m de 0.05 a 0.15 m: 8 J", lambda: trabajo(800 * x, 0.05, 0.15) == 8),
    ("malacate, 400 kg y cable de 1.5 kg/m, 30 m: 124341.75 J",
     lambda: trabajo(_malacate, 0, 30) == sp.Rational("124341.75")),
    ("malacate con 0.05 kWh de entrada: η = 0.6908",
     lambda: cifras(eficiencia(trabajo(_malacate, 0, 30), a_joules(0.05, "kWh")), 4) == "0.6908"),
    ("bombeo de un tanque de radio 1 m con 3 m de agua hasta el borde: 138685.6 J",
     lambda: round(float(trabajo(1000 * _g * sp.pi * (3 - x), 0, 3)), 1) == 138685.6),
    ("η > 1 se rechaza (violaría la conservación de la energía)", lambda: _rechaza(lambda: eficiencia(150, 100))),
    ("x0 = x1 se rechaza", lambda: _rechaza(lambda: trabajo(5 * x, 1, 1))),
]
for nombre, prueba in PRUEBAS_7:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 6.7).** Calcula a mano el trabajo para estirar un resorte de $k=500$ N/m desde su longitud natural ($x=0$) hasta $x=0.2$ m y escribe el número en J.

# %%
mi_valor = None     # escribe un número

if mi_valor is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    try:
        mio = a_numero(mi_valor)
    except ValueError as err:
        print("Revisa tu número:", err)
    else:
        ref = trabajo(500 * x, 0, sp.Rational(1, 5))
        print("La calculadora da: W =", muestra(ref), "J")
        print("coinciden a 3 cifras" if coincide(mio, ref) else
              "NO coinciden: W = ∫_0^0.2 500x dx = ½ k x²; revisa que elevaste 0.2 al cuadrado")
