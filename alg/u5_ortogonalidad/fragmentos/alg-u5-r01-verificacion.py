# ID: ALG-U5-R01
# Libro: Alg. U5 · problema resuelto ALG-U5-01 (verificacion)
# Repositorio: alg/u5_ortogonalidad/fragmentos/alg-u5-r01-verificacion.py
import numpy as np

q = np.array([1, 1, 1, 0, 0])                 # consulta: motor, corriente, par
D = {"d1": [4, 3, 2, 0, 0], "d2": [8, 6, 4, 0, 0],
     "d3": [1, 0, 0, 3, 4], "d4": [0, 0, 0, 2, 1]}
for k, d in D.items():
    d = np.array(d)
    dist = np.linalg.norm(q - d)
    cos = q @ d / (np.linalg.norm(q) * np.linalg.norm(d))
    dist_unit = np.linalg.norm(q / np.linalg.norm(q) - d / np.linalg.norm(d))
    print(f"{k}: distancia {dist:.2f}  coseno {cos:.3f}  "
          f"angulo {np.degrees(np.arccos(cos)):.1f}  "
          f"|q^ - d^|^2 = {dist_unit**2:.3f} = 2 - 2cos = {2 - 2*cos:.3f}")
