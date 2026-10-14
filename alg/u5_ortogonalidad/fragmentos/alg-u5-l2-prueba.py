# ID: ALG-U5-L2
# Libro: Alg. U5 · Laboratorio del Error 5.2 (prueba con codigo)
# Repositorio: alg/u5_ortogonalidad/fragmentos/alg-u5-l2-prueba.py
import numpy as np

x = np.arange(1, 7)                                  # cientos de usuarios
y = np.array([2.6, 4.4, 6.5, 8.6, 10.4, 12.5])      # tiempo de respuesta (ms)
coef = np.polyfit(x, y, 1)
print("polyfit devuelve [pendiente, ordenada] =", np.round(coef, 4))
b, m = coef                                          # error: orden al reves
mal = lambda t: b + m * t                            # "y = 1.99 + 0.54 x"
bien = lambda t: np.polyval(coef, t)
for t in (1, 6, 20):
    print(f"x = {t:2d}: dato {y[t-1] if t <= 6 else '--'}, "
          f"modelo invertido {mal(t):6.2f}, modelo correcto {bien(t):6.2f}")
