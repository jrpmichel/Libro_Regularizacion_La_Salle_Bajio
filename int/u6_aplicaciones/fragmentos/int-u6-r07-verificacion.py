# ID: INT-U6-R07
# Libro: Int. U6 · problema resuelto INT-U6-07 (verificación)
# Repositorio: int/u6_aplicaciones/fragmentos/int-u6-r07-verificacion.py
import sympy as sp

y = sp.symbols("y", real=True)
g, m, lam, H = 9.81, 400, 1.5, 30
W_carga = m * g * H
W_cable = float(sp.integrate(lam * g * (H - y), (y, 0, H)))
W, E = W_carga + W_cable, 0.05 * 3.6e6          # J; 0.05 kWh en J
print(f"carga {W_carga:.0f} J, cable {W_cable:.2f} J, total {W:.2f} J")
print(f"eficiencia {W / E:.3f}; sin el cable {W_carga / E:.3f}")
print(f"cable contado como lam g H por H: {lam * g * H * H:.1f} J")
for Hp in (30, 300):
    c = lam * g * Hp**2 / 2
    print(f"H = {Hp} m: cable {100 * c / (m * g * Hp + c):.1f} %")
