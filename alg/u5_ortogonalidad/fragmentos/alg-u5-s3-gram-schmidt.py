# ID: ALG-U5-S3
# Libro: Alg. U5 · subtema 5.3: Gram-Schmidt paso a paso y factorizacion QR
# Repositorio: alg/u5_ortogonalidad/fragmentos/alg-u5-s3-gram-schmidt.py
import numpy as np
import sympy as sp

def gram_schmidt(vs):
    us = []
    for v in vs:
        u = v - sum((v.dot(w) / w.dot(w) * w for w in us), sp.zeros(len(v), 1))
        if u != sp.zeros(len(v), 1):        # u = 0: v dependia de los anteriores
            us.append(u)
    return us

vs = [sp.Matrix(v) for v in ([1, 1, 0], [1, 0, 1], [0, 1, 1])]
us = gram_schmidt(vs)
print("u:", [list(u) for u in us])
U = sp.Matrix.hstack(*us)
print("U^T U =", (U.T * U).tolist())       # diagonal: ortogonales
Q, R = np.linalg.qr(np.array(sp.Matrix.hstack(*vs), dtype=float))
print("|u_j|   :", [round(float(u.norm()), 4) for u in us])
print("diag R  :", np.round(np.diag(R), 4))   # mismos valores; numpy puede
                                              # cambiar signos
