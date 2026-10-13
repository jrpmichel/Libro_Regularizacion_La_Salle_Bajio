# ID: ALG-U4-L2
# Libro: Alg. U4 · Laboratorio del Error 4.2 (prueba con codigo)
# Repositorio: alg/u4_espacios/fragmentos/alg-u4-l2-prueba.py
def I(T):                                   # transmisor: 4 mA a 0 grados C
    return 4 + 0.16 * T                     # y 20 mA a 100 grados C

print("I(0) =", I(0), "mA")                 # 4: no manda 0 a 0
print("superposicion I(30) + I(50) =", round(I(30) + I(50), 2), "mA")
print("valor real    I(80)         =", round(I(80), 2), "mA")
dI = lambda T: I(T) - I(0)                  # parte lineal: desviacion
print("con desviaciones: 4 +", round(dI(30) + dI(50), 2), "=",
      round(4 + dI(30) + dI(50), 2), "mA")
print("doble de temperatura: I(100) =", I(100), " 2 I(50) =", 2 * I(50))
