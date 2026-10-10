# ID: ALG-U1-S3
# Libro: Alg. U1 · subtema 1.3: producto punto, angulo y proyeccion
# Repositorio: alg/u1_vectores/fragmentos/alg-u1-s3-punto.py
import numpy as np

a = np.array([4.0, 1.0, 1.0])
b = np.array([2.0, 2.0, 1.0])
p = np.dot(a, b)                             # 8 + 2 + 1 = 11
cos_t = p / (np.linalg.norm(a) * np.linalg.norm(b))
cos_t = np.clip(cos_t, -1, 1)                # recorta el redondeo
t = np.arccos(cos_t)                         # arccos devuelve radianes
print("a . b =", p)
print("angulo =", round(t, 4), "rad =", round(np.degrees(t), 2), "grados")

comp = p / np.linalg.norm(b)        # componente escalar de a sobre b
par = (p / np.dot(b, b)) * b        # proyeccion: parte de a paralela a b
perp = a - par                      # lo que sobra es perpendicular a b
print("comp_b a =", round(comp, 4))                # 3.6667
print("a paralela =", np.round(par, 4))
print("a perpendicular =", np.round(perp, 4))
print("perp . b = 0:", abs(np.dot(perp, b)) < 1e-12)   # comprobacion
