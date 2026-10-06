# Código del libro de regularización en matemáticas para ingeniería

Notebooks de Python (Google Colab) y código de los ejemplos del *Libro de regularización*: cálculo diferencial, cálculo integral y álgebra lineal para alumnos que ingresan a las ingenierías de la Universidad La Salle Bajío. Autor: Dr. Jorge Ramón Parra Michel.

## Notebooks disponibles

| Bloque | Unidad | Colab | Notebook | Código |
|---|---|---|---|---|
| Preliminares | Primer notebook | [![Abrir en Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jrpmichel/Libro_Regularizacion_La_Salle_Bajio/blob/main/pre/u0_primer_notebook.ipynb) | [`u0_primer_notebook.ipynb`](pre/u0_primer_notebook.ipynb) | [`pre/u0_primer_notebook/`](pre/u0_primer_notebook) |
| Cálculo diferencial | Dif. U1 · Funciones | [![Abrir en Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/jrpmichel/Libro_Regularizacion_La_Salle_Bajio/blob/main/dif/u1_funciones.ipynb) | [`u1_funciones.ipynb`](dif/u1_funciones.ipynb) | [`dif/u1_funciones/`](dif/u1_funciones) |

Las unidades restantes se publican conforme se terminan. Empieza por el primer notebook (Preliminares).

## Cómo usarlos en Colab

1. Pulsa el botón **Abrir en Colab** de la unidad. Necesitas una cuenta de Google; no se instala nada.
2. Ejecuta las celdas en orden: menú *Entorno de ejecución → Ejecutar todas*.
3. En cada sección, haz a mano la celda **Contrasta** antes de ejecutarla y compara con la calculadora.

Para conservar tus cambios, usa *Archivo → Guardar una copia en Drive* (o el botón *Copiar en Drive*). El notebook de este repositorio no cambia cuando tú lo modificas.

## Cómo descargar el código

- **Todo el repositorio:** botón verde **Code → Download ZIP**.
- **Un archivo:** ábrelo en GitHub y pulsa el icono de descarga (*Download raw file*).
- **Un notebook desde Colab:** *Archivo → Descargar → Descargar .ipynb* (o *.py*).

## Qué hay en cada carpeta

```
<bloque>/                      pre (preliminares), dif, int o alg
  u1_funciones.ipynb           notebook de la unidad (se abre en Colab)
  u1_funciones/
    NN_<subtema>.py            código de cada sección del notebook, en orden
    fragmentos/                programas que el libro imprime (autocontenidos)
herramientas/                  scripts para ensamblar y probar los notebooks
```

Cada archivo empieza con una cabecera de tres líneas: `# ID:` (etiqueta única, por ejemplo `DIF-U1-S1`), `# Libro:` (dónde aparece) y `# Repositorio:` (su ruta aquí). El libro cita el nombre del archivo en el título de cada bloque de código.

## Ejecutarlo en tu computadora

Requiere Python 3.10 o superior.

```
pip install -r requirements.txt
jupyter lab dif/u1_funciones.ipynb
```

Si algo no funciona o encuentras un error en el libro o en el código, abre un *issue* en este repositorio.
