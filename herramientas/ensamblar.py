"""Ensambla un notebook .ipynb a partir de los .py en formato percent de una unidad.

Uso (desde cualquier carpeta):
    python herramientas/ensamblar.py dif/u1_funciones              # escribe dif/u1_funciones.ipynb
    python herramientas/ensamblar.py dif/u1_funciones --ejecutar   # además lo ejecuta y falla si hay error
    python herramientas/ensamblar.py --todas                       # todas las unidades que existan

Convención de los .py: "# %%" abre una celda de código; "# %% [markdown]" abre una de texto
(las líneas van comentadas con "# "). Las líneas de cabecera antes de la primera marca se ignoran.
Los archivos NN_*.py se concatenan en orden (00_, 01_, ..., 99_). La carpeta "fragmentos/" no entra.

Si repo.json trae usuario y repositorio reales, la primera celda lleva el botón "Abrir en Colab".
El notebook se entrega sin salidas y sin metadatos de widgets (GitHub no muestra un notebook
cuyos metadatos de widgets no traen estado)."""
import argparse, json, pathlib, sys
import nbformat
from nbformat.v4 import new_notebook, new_code_cell, new_markdown_cell

RAIZ = pathlib.Path(__file__).resolve().parent.parent
COLAB = 'https://colab.research.google.com/github/{u}/{r}/blob/{b}/{ruta}'
BADGE = '[![Abrir en Colab](https://colab.research.google.com/assets/colab-badge.svg)]({url})'


def repo():
    f = RAIZ / 'repo.json'
    return json.loads(f.read_text(encoding='utf-8')) if f.exists() else {}


def configurado(cfg):
    return all(cfg.get(k) and not cfg[k].startswith('TU_') for k in ('usuario', 'repositorio'))


def url_colab(cfg, ruta):
    return COLAB.format(u=cfg['usuario'], r=cfg['repositorio'], b=cfg.get('rama', 'main'), ruta=ruta)


def leer_percent(ruta):
    celdas, tipo, buf = [], None, []
    def cerrar():
        if tipo is None: return
        texto = '\n'.join(buf).strip('\n')
        if tipo == 'markdown':
            texto = '\n'.join(l[2:] if l.startswith('# ') else l.lstrip('#') for l in texto.split('\n'))
            celdas.append(new_markdown_cell(texto))
        else:
            celdas.append(new_code_cell(texto))
    for linea in ruta.read_text(encoding='utf-8').split('\n'):
        if linea.startswith('# %%'):
            cerrar(); buf = []
            tipo = 'markdown' if '[markdown]' in linea else 'code'
        elif tipo is not None:
            buf.append(linea)
    cerrar()
    return celdas


def ensamblar(unidad, ejecutar=False, salida=None):
    carpeta = RAIZ / unidad
    archivos = sorted(carpeta.glob('[0-9][0-9]_*.py'))
    if not archivos: sys.exit(f'No hay archivos NN_*.py en {carpeta}')
    nb = new_notebook()
    cfg = repo()
    if configurado(cfg):
        nb.cells.append(new_markdown_cell(BADGE.format(url=url_colab(cfg, unidad + '.ipynb'))))
    for f in archivos: nb.cells += leer_percent(f)
    nombre = pathlib.Path(unidad).name
    nb.metadata = {'colab': {'name': f'{nombre}.ipynb', 'provenance': [], 'toc_visible': True},
                   'kernelspec': {'display_name': 'Python 3', 'language': 'python', 'name': 'python3'},
                   'language_info': {'name': 'python'}}
    destino = pathlib.Path(salida) if salida else RAIZ / f'{unidad}.ipynb'
    if ejecutar:
        from nbclient import NotebookClient
        NotebookClient(nb, timeout=300, kernel_name='python3').execute()
        print('Ejecución sin errores.')
        for o in nb.cells[-1].get('outputs', []):
            print(o.get('text', '').strip())
    for i, c in enumerate(nb.cells):      # se entrega sin salidas ni marcas de ejecución
        if c.cell_type == 'code': c.outputs, c.execution_count = [], None
        c.metadata.pop('execution', None)
        c.id = f'celda-{i:03d}'           # ids fijos: reensamblar no cambia el archivo si no cambia el contenido
    nb.metadata.pop('widgets', None)      # evita "Invalid Notebook" en GitHub
    nbformat.validate(nb)
    nbformat.write(nb, destino)
    print(f'{len(nb.cells)} celdas -> {destino.relative_to(RAIZ) if destino.is_relative_to(RAIZ) else destino}')


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('unidad', nargs='?', help='por ejemplo dif/u1_funciones')
    ap.add_argument('-o', '--salida')
    ap.add_argument('--ejecutar', action='store_true')
    ap.add_argument('--todas', action='store_true')
    a = ap.parse_args()
    if a.todas:
        for p in sorted(RAIZ.glob('*/u[0-9]*_*/00_*.py')):
            ensamblar(f'{p.parent.parent.name}/{p.parent.name}', a.ejecutar)
    elif a.unidad:
        ensamblar(a.unidad.strip('/\\').replace('\\', '/'), a.ejecutar, a.salida)
    else:
        ap.error('indica una unidad o usa --todas')


if __name__ == '__main__':
    main()
