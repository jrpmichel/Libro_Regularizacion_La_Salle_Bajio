# ID: ALG-U5-R04
# Libro: Alg. U5 · problema resuelto ALG-U5-04 (verificacion)
# Repositorio: alg/u5_ortogonalidad/fragmentos/alg-u5-r04-verificacion.py
import numpy as np

P = np.arange(6.0)                                     # presion (cientos de kPa)
v = np.array([0.502, 1.278, 2.023, 2.718, 3.381, 3.998])  # salida (V)
for grado in (1, 2):
    A = np.vander(P, grado + 1, increasing=True)       # columnas 1, P, P^2
    coef = np.linalg.solve(A.T @ A, A.T @ v)           # ecuaciones normales
    r = v - A @ coef
    control = np.polyval(coef[::-1], 2.5)              # 250 kPa
    print(f"grado {grado}: coef {np.round(coef, 6)}")
    print(f"   residuos (mV) {np.round(1000 * r, 1)}")
    print(f"   residuo maximo {1000 * abs(r).max():.1f} mV")
    print(f"   control 250 kPa: {control:.4f} V (medido 2.376 V)")
span = 5 * np.polyfit(P, v, 1)[0]                      # escala de la recta
print(f"escala {span:.3f} V; 0.5 % = {5 * span:.1f} mV")
c0, c1, c2 = np.polyfit(P, v, 2)[::-1]
print("p para 2.5 V:", np.round(100 * np.roots([c2, c1, c0 - 2.5]), 1), "kPa")
