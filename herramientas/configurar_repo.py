"""Fija el usuario y el nombre del repositorio de GitHub y regenera README y notebooks.

Uso:  python herramientas/configurar_repo.py <usuario> <repositorio> [rama]
Ejemplo:  python herramientas/configurar_repo.py jrpmichel libro-regularizacion-codigo

Los botones "Abrir en Colab" necesitan la dirección completa del repositorio en GitHub, por eso
se configura una vez, después de crearlo. Escribe repo.json, regenera los README y reensambla
todos los notebooks (sin ejecutarlos)."""
import json, pathlib, subprocess, sys

RAIZ = pathlib.Path(__file__).resolve().parent.parent
if len(sys.argv) not in (3, 4): sys.exit(__doc__)
cfg = {'usuario': sys.argv[1], 'repositorio': sys.argv[2], 'rama': sys.argv[3] if len(sys.argv) == 4 else 'main'}
(RAIZ / 'repo.json').write_text(json.dumps(cfg, indent=2, ensure_ascii=False) + '\n', encoding='utf-8')
aqui = pathlib.Path(__file__).resolve().parent
subprocess.run([sys.executable, str(aqui / 'generar_indice.py')], check=True)
subprocess.run([sys.executable, str(aqui / 'ensamblar.py'), '--todas'], check=True)
print('Listo: sube la carpeta completa al repositorio.')
