"""Genera README.md (raíz) y el README.md de cada unidad a partir de lo que hay en las carpetas.

Uso:  python herramientas/generar_indice.py

Descubre las unidades (<bloque>/u<N>_<slug>/00_*.py), toma el título de la primera celda de texto
de 00_*.py y la descripción de cada archivo de la línea 2 de su cabecera. Los enlaces a Colab
usan repo.json (usuario, repositorio, rama)."""
import json, pathlib, re

RAIZ = pathlib.Path(__file__).resolve().parent.parent
BLOQUES = {'pre': 'Preliminares', 'dif': 'Cálculo diferencial', 'int': 'Cálculo integral', 'alg': 'Álgebra lineal'}
BADGE = '[![Abrir en Colab](https://colab.research.google.com/assets/colab-badge.svg)]({url})'


def cfg():
    return json.loads((RAIZ / 'repo.json').read_text(encoding='utf-8'))


def colab(ruta):
    c = cfg()
    return f'https://colab.research.google.com/github/{c["usuario"]}/{c["repositorio"]}/blob/{c.get("rama", "main")}/{ruta}'


def titulo(carpeta):
    for linea in sorted(carpeta.glob('00_*.py'))[0].read_text(encoding='utf-8').split('\n'):
        m = re.match(r'# # (.+)$', linea)
        if m: return m.group(1).strip()
    return carpeta.name


def descripcion(archivo):
    linea = archivo.read_text(encoding='utf-8').split('\n')[1]
    return linea.split('·', 1)[1].strip() if '·' in linea else ''


def unidades():
    for b in BLOQUES:
        for u in sorted((RAIZ / b).glob('u[0-9]*_*')):
            if u.is_dir() and list(u.glob('00_*.py')): yield b, u


def readme_unidad(b, u):
    ruta = f'{b}/{u.name}'
    L = [f'# {titulo(u)}', '',
         f'{BADGE.format(url=colab(ruta + ".ipynb"))}', '',
         f'Notebook completo: [`{u.name}.ipynb`](../{u.name}.ipynb). Cada celda de código del notebook viene de uno de los archivos `NN_*.py` de esta carpeta, en orden. '
         'Los archivos de [`fragmentos/`](fragmentos) son los programas que el libro imprime; cada uno corre solo, sin el notebook.', '',
         '## Archivos del notebook', '', '| Archivo | Contenido |', '|---|---|']
    for f in sorted(u.glob('[0-9][0-9]_*.py')):
        L.append(f'| [`{f.name}`]({f.name}) | {descripcion(f)} |')
    L += ['', '## Fragmentos que imprime el libro', '', '| Archivo | Dónde aparece |', '|---|---|']
    for f in sorted((u / 'fragmentos').glob('*.py')):
        L.append(f'| [`{f.name}`](fragmentos/{f.name}) | {descripcion(f)} |')
    L += ['', 'Para ejecutar un fragmento en tu computadora: `python fragmentos/<archivo>.py`. En Colab, copia su contenido en una celda de código.', '']
    (u / 'README.md').write_text('\n'.join(L), encoding='utf-8')


def readme_raiz():
    L = ['# Código del libro de regularización en matemáticas para ingeniería', '',
         'Notebooks de Python (Google Colab) y código de los ejemplos del *Libro de regularización*: cálculo diferencial, cálculo integral y álgebra lineal para alumnos que ingresan a las ingenierías de la Universidad La Salle Bajío. Autor: Dr. Jorge Ramón Parra Michel.', '',
         '## Notebooks disponibles', '', '| Bloque | Unidad | Colab | Notebook | Código |', '|---|---|---|---|---|']
    for b, u in unidades():
        ruta = f'{b}/{u.name}'
        L.append(f'| {BLOQUES[b]} | {titulo(u)} | {BADGE.format(url=colab(ruta + ".ipynb"))} | [`{u.name}.ipynb`]({ruta}.ipynb) | [`{ruta}/`]({ruta}) |')
    L += ['', 'Las unidades restantes se publican conforme se terminan. Empieza por el primer notebook (Preliminares).', '',
          '## Cómo usarlos en Colab', '',
          '1. Pulsa el botón **Abrir en Colab** de la unidad. Necesitas una cuenta de Google; no se instala nada.',
          '2. Ejecuta las celdas en orden: menú *Entorno de ejecución → Ejecutar todas*.',
          '3. En cada sección, haz a mano la celda **Contrasta** antes de ejecutarla y compara con la calculadora.', '',
          'Para conservar tus cambios, usa *Archivo → Guardar una copia en Drive* (o el botón *Copiar en Drive*). El notebook de este repositorio no cambia cuando tú lo modificas.', '',
          '## Cómo descargar el código', '',
          '- **Todo el repositorio:** botón verde **Code → Download ZIP**.',
          '- **Un archivo:** ábrelo en GitHub y pulsa el icono de descarga (*Download raw file*).',
          '- **Un notebook desde Colab:** *Archivo → Descargar → Descargar .ipynb* (o *.py*).', '',
          '## Qué hay en cada carpeta', '',
          '```',
          '<bloque>/                      pre (preliminares), dif, int o alg',
          '  u1_funciones.ipynb           notebook de la unidad (se abre en Colab)',
          '  u1_funciones/',
          '    NN_<subtema>.py            código de cada sección del notebook, en orden',
          '    fragmentos/                programas que el libro imprime (autocontenidos)',
          'herramientas/                  scripts para ensamblar y probar los notebooks',
          '```', '',
          'Cada archivo empieza con una cabecera de tres líneas: `# ID:` (etiqueta única, por ejemplo `DIF-U1-S1`), `# Libro:` (dónde aparece) y `# Repositorio:` (su ruta aquí). El libro cita el nombre del archivo en el título de cada bloque de código.', '',
          '## Ejecutarlo en tu computadora', '',
          'Requiere Python 3.10 o superior.', '',
          '```',
          'pip install -r requirements.txt',
          'jupyter lab dif/u1_funciones.ipynb',
          '```', '',
          'Si algo no funciona o encuentras un error en el libro o en el código, abre un *issue* en este repositorio.', '']
    (RAIZ / 'README.md').write_text('\n'.join(L), encoding='utf-8')


if __name__ == '__main__':
    for b, u in unidades():
        readme_unidad(b, u)
    readme_raiz()
    print('README.md y READMEs de unidad generados')
