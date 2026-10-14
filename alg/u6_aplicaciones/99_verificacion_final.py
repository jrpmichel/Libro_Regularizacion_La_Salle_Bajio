# ID: ALG-U6-NB99
# Notebook: alg/u6_aplicaciones.ipynb · verificación final
# Repositorio: alg/u6_aplicaciones/99_verificacion_final.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## Verificación final
# Ejecuta esta celda al terminar. Debe imprimir que todas las pruebas pasaron.

# %%
_ESPERADAS = {1: 9, 2: 8, 3: 7, 4: 8, 5: 8, 6: 8}          # pruebas por sección (el libro cita estos números)
fallas, total = [], 0
for k, pruebas in enumerate([PRUEBAS_1, PRUEBAS_2, PRUEBAS_3, PRUEBAS_4, PRUEBAS_5, PRUEBAS_6], 1):
    assert len(pruebas) == _ESPERADAS[k], (f"la sección {k} tiene {len(pruebas)} pruebas y deben ser "
                                           f"{_ESPERADAS[k]}: falta o sobra alguna")
    resultados = evalua(pruebas)                      # cada caso se ejecuta una sola vez
    print(f"Sección {k}: {sum(ok for _, ok in resultados)}/{len(pruebas)}")
    fallas += [f"sección {k}: {nombre}" for nombre, ok in resultados if not ok]
    total += len(pruebas)
for falla in fallas:
    print("FALLA", falla)
assert not fallas, f"Fallan {len(fallas)} de {total} pruebas (están arriba)"
print(f"Todas las verificaciones pasaron: {total}/{total}")
