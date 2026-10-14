# ID: ALG-U5-R02
# Libro: Alg. U5 · problema resuelto ALG-U5-02 (verificacion)
# Repositorio: alg/u5_ortogonalidad/fragmentos/alg-u5-r02-verificacion.py
import numpy as np

th = np.radians(np.arange(0, 360, 60))          # seis muestras por ciclo
i = np.array([11.0, 8.2, 4.2, -9.0, -6.2, -2.2])  # corriente medida (A)
c, s = np.cos(th), np.sin(th)
print("c.s =", round(c @ s, 12), " |c|^2 =", round(c @ c, 12),
      " |s|^2 =", round(s @ s, 12))
a, b = i @ c / (c @ c), i @ s / (s @ s)         # suma de proyecciones
print(f"a = {a:.3f} A, b = {b:.3f} A, amplitud {np.hypot(a, b):.2f} A, "
      f"fase {np.degrees(np.arctan2(b, a)):.1f} grados")
p = a * c + b * s
r = i - p
print("p =", np.round(p, 3), "\nr =", np.round(r, 3))
print("r.c =", round(r @ c, 12), " r.s =", round(r @ s, 12))
alt = (-1.0) ** np.arange(6)
print("CD:", round(r.mean(), 3), " tercera armonica:", round(r @ alt / 6, 3))
print("|i|^2 =", round(i @ i, 2), " |p|^2 + |r|^2 =", round(p @ p + r @ r, 2))
