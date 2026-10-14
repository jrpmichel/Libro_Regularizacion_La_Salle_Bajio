# ID: ALG-U6-R02
# Libro: Alg. U6 · problema resuelto ALG-U6-02 (verificacion)
# Repositorio: alg/u6_aplicaciones/fragmentos/alg-u6-r02-verificacion.py
import numpy as np

R1 = R2 = R3 = 350.0                     # ohm
Rm, Vex, GF, eps = 10e3, 5.0, 2.0, 1000e-6
Rx = 350.0 * (1 + GF * eps)              # galga deformada: 350.7 ohm
Z = np.array([[R1 + R2, -R1, -R2],
              [-R1, R1 + R3 + Rm, -Rm],
              [-R2, -Rm, R2 + Rm + Rx]])
i = np.linalg.solve(Z, [Vex, 0, 0])
print("i (mA):", np.round(i * 1e3, 4))
Vm = Rm * (i[1] - i[2])                  # V_D - V_B en el detector
aprox = Vex * GF * eps / 4               # formula de la hoja de datos
print(f"detector: {Vm * 1e3:.4f} mV, formula: {aprox * 1e3:.4f} mV")
Vab = Vex * (Rx / (R3 + Rx) - R2 / (R1 + R2))   # sin detector (Rm infinita)
Rbd = R1 * R2 / (R1 + R2) + R3 * Rx / (R3 + Rx)  # vista entre B y D
print(f"sin carga: {Vab * 1e3:.4f} mV;  R_BD = {Rbd:.2f} ohm")
print(f"con el divisor: {Vab * Rm / (Rm + Rbd) * 1e3:.4f} mV")
for nombre, v in (("no linealidad", Vab), ("total", Vm)):
    print(nombre, f"{(v / aprox - 1) * 100:.2f} %")
