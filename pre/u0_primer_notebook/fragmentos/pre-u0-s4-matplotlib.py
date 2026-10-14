# ID: PRE-U0-S4
# Libro: Anexo C · matplotlib: grafica con ejes rotulados y estilos de linea
# Repositorio: pre/u0_primer_notebook/fragmentos/pre-u0-s4-matplotlib.py
import numpy as np
import matplotlib.pyplot as plt

t = np.linspace(0, 2, 200)                         # s
fig, ax = plt.subplots(figsize=(6, 3.5))
for tau, estilo in ((0.5, "k-"), (1.0, "k--")):    # estilo, no solo color
    ax.plot(t, 5 * np.exp(-t / tau), estilo, label=f"tau = {tau} s")
ax.set_xlabel("tiempo t (s)")
ax.set_ylabel("voltaje v (V)")
ax.set_title("Descarga de un capacitor de 5 V")
ax.grid(True, alpha=0.3)
ax.legend()
fig.savefig("descarga.pdf")                   # vectorial, para un reporte
plt.show()
