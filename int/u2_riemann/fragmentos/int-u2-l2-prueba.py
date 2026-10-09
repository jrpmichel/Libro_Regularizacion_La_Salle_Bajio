# ID: INT-U2-L2
# Libro: Int. U2 · Laboratorio del Error 2.2 (prueba con código)
# Repositorio: int/u2_riemann/fragmentos/int-u2-l2-prueba.py
def simpson_ia(y, h):               # la función del asistente: no revisa n
    return h / 3 * (y[0] + y[-1] + 4 * sum(y[1:-1:2]) + 2 * sum(y[2:-1:2]))

Q = [40, 42, 45, 44, 41, 38]        # L/min, cada 2 min: 5 franjas
print("función de la IA:", simpson_ia(Q, 2), "L")                     # 396
print("con una constante:", round(simpson_ia([1] * 6, 2), 4), "(debe dar 10)")
trap = 2 * (sum(Q) - (Q[0] + Q[-1]) / 2)
mixta = simpson_ia(Q[:5], 2) + 2 * (Q[4] + Q[5]) / 2    # Simpson + trapecio
print("trapecio:", trap, "L")
print("Simpson en 4 franjas + trapecio:", round(mixta, 2), "L")
