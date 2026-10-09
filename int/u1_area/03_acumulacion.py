# ID: INT-U1-NB03
# Notebook: int/u1_area.ipynb · sección 1.3 área como acumulación
# Repositorio: int/u1_area/03_acumulacion.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## 3. Área como acumulación: de la velocidad a la posición
#
# Con $v(t)$ dada por una tabla y unida por rectas:
#
# * **Desplazamiento** (área con signo): $\Delta s=A_+-A_-$. La posición final es $s_0+\Delta s$.
# * **Distancia recorrida** (área total): $d=A_++A_-$, el área bajo $|v|$. Siempre $d\geq|\Delta s|$.
#
# Si la velocidad cambia de signo entre dos lecturas, el cruce por cero está en $t_c=t_k+\dfrac{v_k\,(t_{k+1}-t_k)}{v_k-v_{k+1}}$ y cada lado se trata por separado. La misma calculadora sirve para cualquier tasa: caudal → volumen, corriente → carga, potencia → energía.

# %%
def integrar_tabla(ts, vs, s0=0.0):
    """Desplazamiento, distancia, posición máxima y mínima y la tabla de posiciones (t, s)
    para una tasa v(t) definida por puntos y unida por rectas."""
    ts, vs = np.asarray(ts, dtype=float), np.asarray(vs, dtype=float)
    if ts.size != vs.size or ts.size < 2:
        raise ValueError("se necesitan al menos dos puntos y el mismo número de tiempos que de valores.")
    if np.any(np.diff(ts) <= 0):
        raise ValueError("los tiempos deben ir de menor a mayor, sin repetirse.")
    puntos = [(ts[0], vs[0])]
    for k in range(ts.size - 1):                       # inserta los cruces por cero
        t0, t1, v0, v1 = ts[k], ts[k + 1], vs[k], vs[k + 1]
        if v0 * v1 < 0:
            puntos.append((t0 + v0 * (t1 - t0) / (v0 - v1), 0.0))
        puntos.append((t1, v1))
    tp, vp = np.array(puntos).T
    trozos = (vp[:-1] + vp[1:]) / 2 * np.diff(tp)     # trapecios con signo; ninguno cruza el eje
    s = s0 + np.concatenate([[0.0], np.cumsum(trozos)])
    return {"desplazamiento": float(trozos.sum()), "distancia": float(np.abs(trozos).sum()),
            "A+": float(trozos[trozos > 0].sum()), "A-": float(-trozos[trozos < 0].sum()),
            "s_max": float(s.max()), "s_min": float(s.min()), "t": tp, "s": s, "v": vp}


def unidad_area(unidad_v, unidad_t):
    """(unidad del área, aviso). Si la tasa es 'algo/tiempo' y el tiempo coincide, el área queda en 'algo'."""
    unidad_v, unidad_t = unidad_v.strip(), unidad_t.strip()
    if "/" in unidad_v:
        num, den = (p.strip() for p in unidad_v.rsplit("/", 1))
        if den == unidad_t:
            return num, None
        return f"{unidad_v}·{unidad_t}", (f"la tasa está en {unidad_v} y el tiempo en {unidad_t}: convierte uno de los dos "
                                           "antes de interpretar el área (por ejemplo, minutos a horas).")
    return "·".join(u for u in (unidad_v, unidad_t) if u), None


def calculadora_velocidad(t_txt, v_txt, s0, unidad_t, unidad_v, n_cifras):
    try:
        r = integrar_tabla(lista_numeros(t_txt, "los tiempos"), lista_numeros(v_txt, "las velocidades"), s0)
    except ValueError as err:
        print("Revisa la entrada:", err); return
    unidad_s, aviso = unidad_area(unidad_v, unidad_t)
    if aviso:
        print("AVISO:", aviso)
    u = f" {unidad_s}" if unidad_s else ""
    print(f"A+ ≈ {cifras(r['A+'], n_cifras)}{u}   A- ≈ {cifras(r['A-'], n_cifras)}{u}")
    print(f"desplazamiento ≈ {cifras(r['desplazamiento'], n_cifras)}{u}   posición final ≈ "
          f"{cifras(s0 + r['desplazamiento'], n_cifras)}{u}")
    print(f"distancia recorrida ≈ {cifras(r['distancia'], n_cifras)}{u}   posición máxima ≈ {cifras(r['s_max'], n_cifras)}{u}"
          f"   ({n_cifras} cifras significativas, redondeado)")
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(9, 3))
    a1.fill_between(r["t"], 0, r["v"], where=r["v"] >= 0, color="0.8")
    a1.fill_between(r["t"], 0, r["v"], where=r["v"] <= 0, hatch="///", facecolor="white", edgecolor="gray")
    a1.plot(r["t"], r["v"], "o-", color="navy"); a1.axhline(0, color="black", lw=0.7)
    a1.set_xlabel(f"t ({unidad_t})"); a1.set_ylabel(f"tasa ({unidad_v})"); a1.set_title("áreas + y -")
    a2.plot(r["t"], r["s"], "o-", color="navy")
    a2.set_xlabel(f"t ({unidad_t})"); a2.set_ylabel(f"acumulado ({unidad_s})"); a2.set_title("acumulado")
    plt.tight_layout(); plt.show()


widgets.interact(calculadora_velocidad,
    t_txt=widgets.Text(value="0, 4", description="t ="),
    v_txt=widgets.Text(value="4, -4", description="v ="),
    s0=widgets.FloatText(value=0, description="s0"),
    unidad_t=widgets.Text(value="s", description="unidad t"),
    unidad_v=widgets.Text(value="m/s", description="unidad v"),
    n_cifras=widgets.IntSlider(value=4, min=1, max=8, description="cifras sig."));

# %%
# Casos de prueba de la sección 3 (resultado conocido)
_agv = integrar_tabla([0, 2, 8, 10, 12, 13, 16, 17], [0, 1.2, 1.2, 0, 0, -0.8, -0.8, 0])
PRUEBAS_3 = [
    ("semilla: v = 2t en [0, 10] da 100 m", lambda: cerca(integrar_tabla([0, 10], [0, 20])["distancia"], 100)),
    ("v = 4 - 2t en [0, 4]: desplazamiento 0, distancia 8", lambda: cerca(integrar_tabla([0, 4], [4, -4])["desplazamiento"], 0)
                                                                  and cerca(integrar_tabla([0, 4], [4, -4])["distancia"], 8)),
    ("tanque: 3 L/s por 60 s y rampa a 0 en 30 s: 225 L", lambda: cerca(integrar_tabla([0, 60, 90], [3, 3, 0])["desplazamiento"], 225)),
    ("AGV: final 6.4, odómetro 12.8, máximo 9.6", lambda: cerca(_agv["desplazamiento"], 6.4) and cerca(_agv["distancia"], 12.8)
                                                       and cerca(_agv["s_max"], 9.6)),
    ("tiempos que no crecen se rechazan", lambda: _rechaza(lambda: integrar_tabla([0, 2, 2], [1, 1, 1]))),
    ("unidades: m/s por s da m; km/h por min avisa", lambda: unidad_area("m/s", "s") == ("m", None)
                                                         and unidad_area("km/h", "min")[1] is not None),
]
for nombre, prueba in PRUEBAS_3:
    print("OK  " if prueba() else "FALLA", nombre)

# %% [markdown]
# **Contrasta (sección 3).** La velocidad de un carrito pasa por $(0,2)$, $(3,2)$, $(4,-2)$ y $(6,-2)$ ($t$ en s, $v$ en m/s), unidos por rectas. Calcula a mano el desplazamiento y la distancia recorrida.

# %%
mi_desplazamiento = None    # escribe un número
mi_distancia = None         # escribe un número

if mi_desplazamiento is None or mi_distancia is None:     # la referencia se muestra después de tu cálculo
    print("falta tu cálculo: escríbelo arriba y vuelve a ejecutar la celda")
else:
    r = integrar_tabla([0, 3, 4, 6], [2, 2, -2, -2])
    print("La calculadora da: desplazamiento =", cifras(r["desplazamiento"], 4), "m   distancia =", cifras(r["distancia"], 4), "m")
    print("coinciden" if cerca(mi_desplazamiento, r["desplazamiento"], 1e-6) and cerca(mi_distancia, r["distancia"], 1e-6) else
          "NO coinciden: la velocidad cruza el cero en t = 3.5 s; separa los dos triángulos que se forman ahí")
