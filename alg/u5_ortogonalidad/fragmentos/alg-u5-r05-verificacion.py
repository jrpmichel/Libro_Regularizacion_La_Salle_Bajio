# ID: ALG-U5-R05
# Libro: Alg. U5 · problema resuelto ALG-U5-05 (verificacion)
# Repositorio: alg/u5_ortogonalidad/fragmentos/alg-u5-r05-verificacion.py
import numpy as np

x, y = np.meshgrid([5, 15, 25, 35], [5, 15, 25])   # centros de celda (m)
z = np.array([[100.50, 100.65, 100.88, 101.34],
              [100.92, 101.35, 101.58, 101.70],
              [101.32, 101.66, 102.06, 102.41]])   # cotas (m)
x, y, z = x.ravel(), y.ravel(), z.ravel()
u, w = x - x.mean(), y - y.mean()                  # columnas centradas
print("1.u, 1.w, u.w:", u.sum(), w.sum(), u @ w)   # 0 0 0: ortogonales
b1, b2 = u @ z / (u @ u), w @ z / (w @ w)
print(f"plano: z = {z.mean():.3f} + {b1:.4f}(x-20) + {b2:.4f}(y-15)")
for baja in (0.0, 0.12):                           # plataforma con 1 % en x
    h = np.mean(z - 0.01 * x) - baja
    d = z - (h + 0.01 * x)                         # + corte, - relleno
    print(f"h = {h:.3f} m: corte {100 * d[d > 0].sum():.0f} m3, "
          f"relleno {-100 * d[d < 0].sum():.0f} m3")
    print(f"   relleno max {-d.min():.2f} m, corte max {d.max():.2f} m")
print(np.round((z - (np.mean(z - 0.01 * x) + 0.01 * x)).reshape(3, 4), 3))
