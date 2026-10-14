# ID: ALG-U6-L3
# Libro: Alg. U6 · Laboratorio del Error: la inversa que solo probo traslaciones
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-l3-prueba.py
import numpy as np

def H(R, t):
    T = np.eye(4); T[:3, :3] = R; T[:3, 3] = t; return T

def inv_ia(T):                           # la funcion de la IA
    return H(T[:3, :3].T, -T[:3, 3])

def inv_bien(T):
    R, t = T[:3, :3], T[:3, 3]
    return H(R.T, -R.T @ t)

Tt = H(np.eye(3), [1, 2, 3])                      # la prueba de la IA
print("prueba de la IA:", np.allclose(Tt @ inv_ia(Tt), np.eye(4)))
R = np.column_stack([[0, -1, 0], [0, 0, -1], [1, 0, 0]])
T = H(R, [-0.3, 0.0, 0.5])                        # camara de ALG-U6-07
print("T inv_ia(T) =\n", np.round(T @ inv_ia(T), 3))
print("base en la camara, IA:", inv_ia(T)[:3, 3] + 0.0,
      " bien:", inv_bien(T)[:3, 3] + 0.0)
print("T inv_bien(T) = I:", np.allclose(T @ inv_bien(T), np.eye(4)))
