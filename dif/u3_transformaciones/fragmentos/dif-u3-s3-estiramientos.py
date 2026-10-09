# ID: DIF-U3-S3
# Libro: Dif. U3 · subtema 3.3: amplitud, periodo y factor de escala
# Repositorio: dif/u3_transformaciones/fragmentos/dif-u3-s3-estiramientos.py
import numpy as np

for b in (0.5, 1, 2):
    T = 2 * np.pi / abs(b)                       # periodo de sen(bx)
    print(f"sen({b}x): periodo = {T:.4f} = {T / np.pi:.2f} pi")
for a in (0.5, 2, -3):
    print(f"{a}·sen x: amplitud {abs(a)}, rango [{-abs(a)}, {abs(a)}]")

x = np.array([1.0, 4.0, 9.0])
print(np.sqrt(x / 4), np.sqrt(x) / 2)           # sqrt(x/4) = sqrt(x)/2
