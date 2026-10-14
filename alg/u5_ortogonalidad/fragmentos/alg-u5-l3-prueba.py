# ID: ALG-U5-L3
# Libro: Alg. U5 · Laboratorio del Error 5.3 (prueba con codigo)
# Repositorio: alg/u5_ortogonalidad/fragmentos/alg-u5-l3-prueba.py
import numpy as np

t = np.arange(1, 7)                                  # mes
d = np.array([120, 135, 128, 150, 158, 149])        # demanda (piezas)
p5 = np.polyfit(t, d, 5)                             # 6 coeficientes, 6 datos
print("residuos grado 5:", np.round(d - np.polyval(p5, t), 6))   # ceros
print("pronostico mes 7, grado 5:", round(np.polyval(p5, 7), 1))
p1 = np.polyfit(t, d, 1)
r1 = d - np.polyval(p1, t)
print("recta:", np.round(p1[::-1], 3), " RMSE:", round(np.sqrt(r1 @ r1 / 6), 1))
print("pronostico mes 7, recta:", round(np.polyval(p1, 7), 1))
p5b = np.polyfit(t[:5], d[:5], 4)                    # validacion: sin el mes 6
p1b = np.polyfit(t[:5], d[:5], 1)
print("mes 6 sin usarlo: grado 4 da", round(np.polyval(p5b, 6), 1),
      "; recta da", round(np.polyval(p1b, 6), 1), "; real", d[5])
