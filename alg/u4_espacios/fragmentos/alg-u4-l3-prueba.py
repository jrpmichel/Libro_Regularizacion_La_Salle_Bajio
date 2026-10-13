# ID: ALG-U4-L3
# Libro: Alg. U4 · Laboratorio del Error 4.3 (prueba con codigo)
# Repositorio: alg/u4_espacios/fragmentos/alg-u4-l3-prueba.py
import numpy as np

def R(grados):
    t = np.deg2rad(grados)
    return np.array([[np.cos(t), -np.sin(t)], [np.sin(t), np.cos(t)]])

pc = np.array([10, 0])                      # en ejes de la camara (cm)
for nombre, M in (("IA, R(-30)", R(-30)), ("correcto, R(30)", R(30))):
    p = M @ pc
    ang = np.degrees(np.arctan2(p[1], p[0]))
    print(f"{nombre}: p = {np.round(p, 2)} cm, direccion {ang:.0f} grados")
print("eje x de la camara en la base:", np.round(R(30) @ [1, 0], 4))
print("de la base a la camara:", np.round(R(-30) @ R(30) @ pc, 2))
