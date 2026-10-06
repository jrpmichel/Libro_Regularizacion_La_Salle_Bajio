# ID: PRE-U0-NB00
# Notebook: pre/u0_primer_notebook.ipynb · primer notebook
# Repositorio: pre/u0_primer_notebook/00_primer_notebook.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# # Primer notebook
# **Libro de regularización en matemáticas para ingeniería · ULSA Bajío**
#
# Este notebook acompaña el capítulo *Python y Google Colab*. No necesitas instalar nada ni saber programar: solo ejecutar las celdas y mirar qué sale.
#
# Hay dos tipos de celdas. Esta es de **texto**. Las que tienen un botón ▶ a la izquierda son de **código**.
#
# **Cómo ejecutarlas.** Pulsa el botón ▶ de una celda de código, o selecciónala y pulsa Mayús+Enter. El resultado aparece justo debajo. Ejecuta las celdas de arriba hacia abajo: cada una usa lo que definieron las anteriores. Para ejecutar todo de una vez: menú *Entorno de ejecución → Ejecutar todas*.

# %% [markdown]
# ## Tabla y gráfica de una función
# La función $f(x)=x^2+1$ asigna a cada $x$ el valor $x^2+1$. Esta celda calcula $f$ para $x=0,1,2,3,4$, imprime la tabla y dibuja los puntos. Todo lo que sigue a un `#` es un comentario: Python lo ignora.

# %%
import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return x**2 + 1          # la regla; cámbiala y vuelve a ejecutar

x = np.arange(0, 5)          # entradas: 0, 1, 2, 3, 4
y = f(x)                     # salidas

print("x =", x)
print("y =", y)

plt.plot(x, y, "o-")
plt.grid(True)
plt.xlabel("x")
plt.ylabel("y")
plt.show()

# %% [markdown]
# ## Pruébalo tú
# Antes de ejecutar cada cambio, **predice** la tabla: calcula a mano $f(2)$ y anótalo. Después ejecuta y compara.
#
# 1. En la celda anterior, cambia `return x**2 + 1` por `return 2*x + 1` y ejecútala otra vez. ¿Los puntos quedan sobre una recta?
# 2. Vuelve a `x**2 + 1` y cambia `np.arange(0, 5)` por `np.arange(-3, 4)`. ¿Qué pasa con la gráfica para $x$ negativos?
#
# Si aparece un mensaje de error en rojo, lee su última línea: casi siempre dice qué se escribió mal. Si dice `NameError`, ejecutaste una celda sin haber ejecutado antes las de arriba.

# %% [markdown]
# ## Comprobación
# Ejecuta esta celda con la regla original ($x^2+1$). Debe imprimir que todo funciona.

# %%
esperado = np.array([1, 2, 5, 10, 17])
assert np.array_equal(f(np.arange(0, 5)), esperado), "La regla f ya no es x**2 + 1: restáurala para esta comprobación."
print("Todo funciona: numpy y matplotlib están listos.")
