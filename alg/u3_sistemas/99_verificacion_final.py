# ID: ALG-U3-NB99
# Notebook: alg/u3_sistemas.ipynb · verificación final
# Repositorio: alg/u3_sistemas/99_verificacion_final.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## Verificación final
# Ejecuta esta celda al terminar. Debe imprimir que todas las pruebas pasaron.

# %%
fallas, total = [], 0
for k, pruebas in enumerate([PRUEBAS_1, PRUEBAS_2, PRUEBAS_3, PRUEBAS_4, PRUEBAS_5], 1):
    resultados = evalua(pruebas)                      # cada caso se ejecuta una sola vez
    print(f"Sección {k}: {sum(ok for _, ok in resultados)}/{len(pruebas)}")
    fallas += [f"sección {k}: {nombre}" for nombre, ok in resultados if not ok]
    total += len(pruebas)
for falla in fallas:
    print("FALLA", falla)
assert not fallas, f"Fallan {len(fallas)} de {total} pruebas (están arriba)"
print(f"Todas las verificaciones pasaron: {total}/{total}")
