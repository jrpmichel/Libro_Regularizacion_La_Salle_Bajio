# Herramientas del repositorio

Scripts para quien mantiene el código. Los alumnos no los necesitan.

| Script | Para qué |
|---|---|
| `ensamblar.py` | Arma `<bloque>/<unidad>.ipynb` a partir de los `NN_*.py` de la unidad. Con `--ejecutar` lo corre completo y falla si una celda da error. Entrega el notebook sin salidas y sin metadatos de widgets. |
| `probar_fragmentos.py` | Audita las cabeceras (`# ID:`, `# Libro:`, `# Repositorio:`) y ejecuta cada fragmento. Dentro de la carpeta del libro compara además con los `\fragmento{...}` del `.tex`. |
| `generar_indice.py` | Genera el `README.md` de la raíz y el de cada unidad. |
| `configurar_repo.py` | Fija usuario y nombre del repositorio de GitHub en `repo.json` y regenera README y notebooks, con los botones de Colab. |

Uso, desde la raíz del repositorio:

```
python herramientas/ensamblar.py dif/u1_funciones --ejecutar
python herramientas/probar_fragmentos.py dif/u1_funciones
python herramientas/configurar_repo.py <usuario> <repositorio>
```

Dependencias: `pip install -r requirements.txt`.

## Etiquetas

La línea 1 de cada archivo es `# ID: ...` en mayúsculas. El nombre del archivo empieza con la misma etiqueta en minúsculas.

| Etiqueta | Qué es | Archivo |
|---|---|---|
| `DIF-U1-NB##` | sección del notebook | `NN_<slug>.py` |
| `DIF-U1-S#` | código del subtema # | `fragmentos/dif-u1-s#-<slug>.py` |
| `DIF-U1-R##` | verificación del resuelto `DIF-U1-##` | `fragmentos/dif-u1-r##-verificacion.py` |
| `DIF-U1-L#` | prueba con código del Laboratorio del Error # | `fragmentos/dif-u1-l#-prueba.py` |

Cabecera de cada fragmento (3 líneas; el libro no las imprime):

```
# ID: DIF-U1-S1
# Libro: Dif. U1 · subtema 1.1: tabla, diferencias y gráfica
# Repositorio: dif/u1_funciones/fragmentos/dif-u1-s1-tabla-a-curva.py
```

## Cómo entra el código al libro

El `.tex` no copia el código: lo lee del archivo con `\fragmento{tipo}{ID}{archivo}`. Cada unidad declara su carpeta con `\renewcommand{\codigodir}{notebooks/dif/u1_funciones/fragmentos}`. El título de la caja imprime `<archivo>.py`, así que el lector ubica el archivo en el repositorio. Para cambiar código impreso, se edita el fragmento y se recompila.

## Formato de los `.py`

Formato percent de jupytext: `# %%` abre una celda de código y `# %% [markdown]` una de texto, con cada línea comentada como `# `. Las líneas anteriores a la primera marca (cabecera) se ignoran. Las líneas con `%pip` solo funcionan en IPython o Colab.

## Reglas del código

- Ejecutable en Colab sin instalar nada más allá de la primera celda.
- Comentarios técnicos y breves; los textos para el alumno van en español de México, con punto decimal.
- Las calculadoras validan la entrada, redondean (no truncan) a las cifras significativas pedidas, muestran ambas raíces con signo y van acompañadas de una celda Contrasta y casos de prueba.
- Cada calculadora nueva suma su caso a la lista `PRUEBAS_n` de la sección.
