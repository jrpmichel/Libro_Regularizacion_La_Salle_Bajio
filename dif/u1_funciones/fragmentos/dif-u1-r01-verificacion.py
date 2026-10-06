# ID: DIF-U1-R01
# Libro: Dif. U1 · problema resuelto DIF-U1-01 (verificación)
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-r01-verificacion.py
def es_funcion(pares):
    vistos = {}
    for t, T in pares:
        if vistos.setdefault(t, T) != T:      # mismo t con otra T
            return False
    return True

A = [(0, 20), (5, 55), (10, 88), (15, 115)]
B = [(0, 20), (5, 55), (5, 60), (10, 88)]
C = [(0, 20), (5, 20), (10, 88), (15, 88)]
print(es_funcion(A), es_funcion(B), es_funcion(C))   # True False True
