# ID: ALG-U1-S1
# Libro: Alg. U1 · subtema 1.1: componentes del vector PQ, magnitud y angulo de direccion
# Repositorio: alg/u1_vectores/fragmentos/alg-u1-s1-componentes.py
import numpy as np
import matplotlib.pyplot as plt

P = np.array([1.0, 1.0])
Q = np.array([4.0, 5.0])
v = Q - P                              # componentes: llegada menos salida
print("PQ =", v)                       # [3. 4.]
print("magnitud =", np.linalg.norm(v))                 # 5.0
ang = np.degrees(np.arctan2(v[1], v[0]))   # arctan2 respeta el cuadrante
print("angulo con el eje x =", round(ang, 2), "grados")   # 53.13

w = np.array([-3.0, -4.0])             # misma magnitud, sentido opuesto
mal = np.degrees(np.arctan(w[1] / w[0]))
bien = np.degrees(np.arctan2(w[1], w[0]))
print("arctan(-4/-3) =", round(mal, 2))        # 53.13: mal
print("arctan2(-4, -3) =", round(bien, 2))     # -126.87

plt.quiver(*P, *v, angles="xy", scale_units="xy", scale=1, color="navy")
plt.plot([P[0], Q[0]], [P[1], P[1]], "k:")     # componente x
plt.plot([Q[0], Q[0]], [P[1], Q[1]], "k:")     # componente y
plt.axis("equal"); plt.xlim(0, 6); plt.ylim(0, 6)
plt.grid(True); plt.xlabel("x"); plt.ylabel("y"); plt.show()
