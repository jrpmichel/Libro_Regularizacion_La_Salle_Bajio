# ID: ALG-U1-NB99
# Notebook: alg/u1_vectores.ipynb · verificación final
# Repositorio: alg/u1_vectores/99_verificacion_final.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## Verificación final
# Ejecuta esta celda al terminar. Debe imprimir que todas las pruebas pasaron.

# %%
todas = PRUEBAS_1 + PRUEBAS_2 + PRUEBAS_3 + PRUEBAS_4 + PRUEBAS_5
fallas = [nombre for nombre, prueba in todas if not prueba()]
assert not fallas, f"Fallan {len(fallas)} pruebas: {fallas}"
print(f"Todas las verificaciones pasaron: {len(todas)}/{len(todas)}")
