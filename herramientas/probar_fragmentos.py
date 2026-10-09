"""Audita y ejecuta el código de una unidad.

Uso:  python herramientas/probar_fragmentos.py dif/u1_funciones

1. Etiquetas. Cada archivo de fragmentos/ lleva tres líneas de cabecera
       # ID: DIF-U1-S1
       # Libro: Dif. U1 · subtema 1.1: <descripción>
       # Repositorio: dif/u1_funciones/fragmentos/dif-u1-s1-<slug>.py
   El ID en minúsculas debe ser el prefijo del nombre del archivo y la ruta del repositorio debe coincidir
   con la ubicación real. Los archivos NN_*.py llevan "# ID: DIF-U1-NBnn".
   Si el repositorio vive dentro de la carpeta del libro (existe ../Secciones/<bloque>-<unidad>.tex),
   cada \\fragmento{tipo}{ID}{archivo} del .tex debe apuntar a un archivo con ese ID y ningún
   fragmento puede quedar sin usar. En un clon del repositorio solo se auditan las cabeceras.
   Los fragmentos solo pueden traer caracteres ASCII o Latin-1 (acentos, ñ, °, ², ±): con XeLaTeX, listings
   imprime fuera de orden los demás (√, π, ≈...).
2. Ejecución. Cada fragmento corre como script independiente (matplotlib sin ventana).
Termina con código 1 si falla cualquier comprobación."""
import os, pathlib, re, subprocess, sys
os.environ['MPLBACKEND'] = 'Agg'

if len(sys.argv) != 2: sys.exit(__doc__)
ruta_unidad = sys.argv[1].strip('/\\').replace('\\', '/')
bloque, unidad = ruta_unidad.split('/')
raiz = pathlib.Path(__file__).resolve().parent.parent      # raíz del repositorio (carpeta notebooks)
carpeta = raiz / bloque / unidad
frags = sorted((carpeta / 'fragmentos').glob('*.py'))
tex = raiz.parent / 'Secciones' / f'{bloque}-{unidad.replace("_", "-")}.tex'
errores = []

# --- 1. etiquetas -----------------------------------------------------------
ids = {}
for f in frags:
    cab = f.read_text(encoding='utf-8').split('\n')[:3]
    m = re.match(r'# ID: ([A-Z]+-U\d+-[A-Z]\d+)$', cab[0])
    if not m: errores.append(f'{f.name}: falta "# ID: ..." en la línea 1'); continue
    tag = m.group(1); ids[tag] = f.stem
    if not f.stem.startswith(tag.lower() + '-') and f.stem != tag.lower():
        errores.append(f'{f.name}: el ID {tag} no es prefijo del nombre')
    if not cab[1].startswith('# Libro: '): errores.append(f'{f.name}: falta "# Libro:" en la línea 2')
    esperado = f'# Repositorio: {bloque}/{unidad}/fragmentos/{f.name}'
    if cab[2] != esperado: errores.append(f'{f.name}: línea 3 debe ser "{esperado}"')
# Fuera de Latin-1 (√, π, ≈...), listings imprime los caracteres fuera de orden con XeLaTeX
for f in frags:
    for n, linea in enumerate(f.read_text(encoding='utf-8').split('\n')[3:], 4):
        raros = sorted({c for c in linea if ord(c) > 0xFF})
        if raros: errores.append(f'{f.name}, línea {n}: {raros} no se imprime bien en el libro; escribe sqrt, pi, aprox.')
for f in sorted(carpeta.glob('[0-9][0-9]_*.py')):
    primera = f.read_text(encoding='utf-8').split('\n')
    if not re.match(rf'# ID: {bloque.upper()}-U\d+-NB{f.name[:2]}$', primera[0]):
        errores.append(f'{f.name}: ID de notebook ausente o incorrecto')
    if primera[2] != f'# Repositorio: {bloque}/{unidad}/{f.name}':
        errores.append(f'{f.name}: línea 3 debe ser "# Repositorio: {bloque}/{unidad}/{f.name}"')
if tex.exists():
    usados = re.findall(r'\\fragmento\{(\w+)\}\{([A-Z0-9-]+)\}\{([\w-]+)\}', tex.read_text(encoding='utf-8'))
    for tipo, tag, arch in usados:
        if ids.get(tag) != arch: errores.append(f'.tex: {tag} apunta a "{arch}" y el archivo con ese ID es "{ids.get(tag)}"')
    for tag, arch in ids.items():
        if arch not in [a for _, _, a in usados]: errores.append(f'{arch}.py no se usa en {tex.name}')
    print(f'{len(usados)} \\fragmento en {tex.name}, {len(frags)} archivos en fragmentos/')
else:
    print(f'Sin {tex.name} junto al repositorio: se omite la comparación con el libro ({len(frags)} archivos en fragmentos/)')

# --- 2. ejecución -----------------------------------------------------------
fallas = 0
for f in frags:
    r = subprocess.run([sys.executable, str(f)], capture_output=True, text=True, timeout=120)
    fallas += r.returncode != 0
    print(f'{"ok" if r.returncode == 0 else "FALLA":6s} {f.name}')
    if r.returncode: print(r.stderr.strip().splitlines()[-1])
for e in errores: print('ETIQUETA:', e)
print(f'{len(frags) - fallas}/{len(frags)} fragmentos sin error; {len(errores)} problemas de etiquetas')
sys.exit(1 if (fallas or errores) else 0)
