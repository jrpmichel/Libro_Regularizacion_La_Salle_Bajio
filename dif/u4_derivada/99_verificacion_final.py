# ID: DIF-U4-NB99
# Notebook: dif/u4_derivada.ipynb · verificación final
# Repositorio: dif/u4_derivada/99_verificacion_final.py
# Formato percent (jupytext): "# %%" abre celda de código, "# %% [markdown]" abre celda de texto.

# %% [markdown]
# ## Verificación final
# Ejecuta esta celda al terminar. Debe imprimir que todas las pruebas pasaron.

# %%
todas = PRUEBAS_1 + PRUEBAS_2 + PRUEBAS_3 + PRUEBAS_4 + PRUEBAS_5 + PRUEBAS_6
fallas = [nombre for nombre, prueba in todas if not prueba()]
assert not fallas, f"Fallan {len(fallas)} pruebas: {fallas}"
print(f"Todas las verificaciones pasaron: {len(todas)}/{len(todas)}")
