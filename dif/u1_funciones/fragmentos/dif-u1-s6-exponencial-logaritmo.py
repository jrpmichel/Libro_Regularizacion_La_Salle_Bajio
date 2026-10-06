# ID: DIF-U1-S6
# Libro: Dif. U1 · subtema 1.6: exponencial y logaritmo
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-s6-exponencial-logaritmo.py
import numpy as np

x = np.array([-2, -1, 0, 1, 2, 3])
print("2^x      :", 2.0**x)              # crece y siempre es positiva
print("(1/2)^x  :", 0.5**x)              # decrece y siempre es positiva
print("ln(e^x)  :", np.log(np.exp(x)))   # ln deshace a exp: devuelve x
y = np.array([0.5, 1, 2, 10])
print("e^(ln y) :", np.exp(np.log(y)))   # exp deshace a ln, y solo para y > 0
x_sol = np.log(7) / np.log(5)            # 5^x = 7  ->  x = ln 7 / ln 5
print("x =", round(x_sol, 4), "| 5^x =", round(5.0**x_sol, 4))   # 1.2091 | 7.0
