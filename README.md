# DNAinvent 🧬

Herramientas para trabajar con secuencias de ADN y ARN: doble hélice, transcripción,
traducción, codificación binaria, alineamiento, permutaciones y edición con CRISPR/Cas9.
Y un **microscopio 3D** para viajar dentro de la secuencia: de la doble hélice a un átomo,
a su núcleo, a un protón y a sus quarks.
Funciona en consola, en catalán, castellano e inglés.

Este repositorio es el **núcleo** (la lógica, la consola y el microscopio de escritorio). La
**web** está en su propio repositorio, `dnainvent_web`, que usa este como paquete de Python.

## Instalación

El microscopio usa la física atómica de **AtomicWorld**, que se instala desde su carpeta.
Con las carpetas así:

```
Documents/
├── project-DNAinvent/DNAinvent
└── project-AtomicWorld/AtomicWorld
```

```bash
pip install -r requirements.txt     # instala también ../../project-AtomicWorld/AtomicWorld
```

Si solo quieres la consola, sin microscopio, basta con `pip install colorama`.

## Uso

```bash
python main.py                 # consola (la primera vez pregunta el idioma y lo guarda)
python main.py --lang ca       # consola en catalán (también: es, en)
python -m pytest               # tests
```

**Idiomas:** la primera vez la consola pregunta el idioma y lo guarda en
`~/.dnainvent/settings.json`; se cambia con la opción **L** del menú o con la tecla **L**
dentro del microscopio. La consola y el microscopio comparten ese archivo.

## Microscopio 3D

Viaje en primera persona dentro de la secuencia, desde la doble hélice hasta los quarks.
Se abre desde la consola (**7] Microscopio 3D**, o al terminar la opción **1] Doble hélice**)
o directamente:

```bash
python -m microscope ATGCGTACGTTAGC
python -m microscope ATGC --lang ca --realistic
python -m microscope <secuencia larga> --start 120 --count 60     # un tramo
```

El microscopio enseña como mucho **60 pares** a la vez (`dnainvent/scene/limits.py`; en la
web, 100). Con una secuencia más larga no se rechaza: se elige desde qué base empezar, y el
HUD sigue numerando los pares como en la secuencia completa.

| Nivel | Escala | Qué se ve |
|---|---|---|
| 1 · Hélice | nm | la molécula completa como suma de nubes electrónicas; puentes de H en amarillo |
| 2 · Átomo | pm | orbitales del átomo y de sus vecinos enlazados; el núcleo a tamaño real |
| 3 · Núcleo | fm | protones y neutrones (didáctico) o densidad nuclear (realista) |
| 4 · Nucleón | 0,1 fm | quarks de valencia, tubo de gluones en Y y mar de partones |

| Tecla | Acción |
|---|---|
| Clic / Esc | capturar / soltar el ratón (Esc otra vez: salir) |
| W A S D · Espacio / Ctrl | moverse · subir / bajar |
| Q / E · Shift · Rueda | girar · turbo ×5 · multiplicador de velocidad |
| **Enter** | piloto automático: bajar un nivel (al objeto de la mira) |
| **Mayús+Enter** | viaje completo automático hasta un quark |
| **Retroceso** | piloto automático: subir un nivel |
| Tab | modo didáctico ↔ realista |
| **L** | idioma: català → English → español (se recuerda) |
| H · F11 | ayuda · pantalla completa |
| **C** | (comparando proteínas) pantalla partida ↔ superpuestas |
| **V** | (pantalla partida) cambiar el lado que manejas |

**Física:**

- **Geometría:** ADN forma B construido par a par desde plantillas atómicas reales (PyMOL,
  con hidrógenos y el fosfato cargado), con el mismo paso que `dnainvent/helix3d.py`:
  **10,5 pares por vuelta** (34,3°) y 3,38 Å por par. Los puentes de hidrógeno no se
  escriben a mano: salen de la geometría (2 en A–T, 3 en G–C) y los tests lo comprueban.
- **Del átomo hacia abajo**, la física es la de **AtomicWorld** (`atomicworld.scene`):
  orbitales hidrogenoides con la carga efectiva de Slater, lóbulos p reales con la regla de
  Hund, núcleo de radio 1,2 fm·A^(1/3) con los nucleones empaquetados (o perfil de
  Woods–Saxon en modo realista) y quarks con las masas del PDG. Los datos de H, C, N, O y P
  (radio, color, isótopo más abundante) vienen de su tabla periódica. Así las dos apps
  enseñan exactamente lo mismo, y una mejora en AtomicWorld llega también aquí.

**Limitaciones:** los orbitales son de átomo aislado (en la molécula se deforman al
enlazarse); la hélice es ideal (el ADN real varía con la secuencia), y al pasar las plantillas
de 10 a 10,5 pares por vuelta el enlace O3'–P entre pares queda en 1,57 Å en lugar de 1,60 Å;
el núcleo es una "gota líquida" sin capas; el movimiento de los quarks y el tubo en Y son una
ilustración: la QCD no tiene trayectorias.

`tools/microscope_render_test.py` hace el viaje completo con el piloto automático, sin
ventana, y guarda capturas PNG de cada nivel.

**Gráfica NVIDIA en Windows:** Configuración → Sistema → Pantalla → Gráficos → añade el
`python.exe` de tu entorno → Alto rendimiento. Al arrancar, la consola dice qué GPU se usa.

## Estructura

```
dnainvent/             lógica biológica (sin input ni print)
  scene/               física del microscopio, sin OpenGL (necesita NumPy y AtomicWorld)
    molecule.py        secuencia (o un tramo) -> todos los átomos de la doble hélice
    elements.py        H, C, N, O y P con los datos de AtomicWorld
    describe.py        resumen de la molécula para la web (/api/structure)
    limits.py          60 pares en escritorio, 100 en la web; tramos (window)
    protein.py         estructura de proteína -> todos los átomos, con la forma de molecule.py
    templates/         plantilla atómica de cada par de bases (tools/generate_templates.py)
  protein/             del ADN a la proteína 3D (en desarrollo, ver más abajo)
  database/            base de datos de genes, variantes y enfermedades reales (ver más abajo)
microscope/            microscopio de escritorio (moderngl, pyglet, pillow): solo dibuja
  __main__.py          python -m microscope ATGC --lang ca
  app.py · engine.py   ventana y teclas · pila de niveles y transiciones
  autopilot.py         Enter / Mayús+Enter / Retroceso
  ship.py · hud.py     nave con velocidad logarítmica · HUD
  levels/              helix, atom, nucleus, nucleon
  render/              puntos, líneas y shaders
i18n/                  traducciones: un JSON por idioma (ca, en, es), compartido por la consola,
                       el microscopio y la web
ui/                    consola (console.py), idioma guardado (settings.py) y lanzador del
                       microscopio en otro proceso (launcher.py)
tests/                 tests automáticos (los del microscopio no necesitan tarjeta gráfica)
tools/                 generate_templates.py y generate_amino_acid_templates.py (PyMOL),
                       microscope_render_test.py y protein_info.py
pyproject.toml         permite instalar dnainvent e i18n como paquete (lo usa dnainvent_web)
```

Para añadir un idioma, copia `i18n/locales/es.json` con otro nombre y traduce los textos: los
tests comprueban que no falte ninguna clave. Los textos de la web (claves `web.*`) y del
microscopio (`microscope.*`) también están aquí, así consola, microscopio y web comparten las
traducciones. Los nombres de los elementos vienen de AtomicWorld (`atomicworld.names`).

## Proteínas 3D (en desarrollo)

De la secuencia de ADN a la proteína plegada: se busca el marco de lectura, se identifica la
proteína (BLAST) y se muestra su estructura de AlphaFold DB, o se predice si es inventada.

| Fase | Qué hay | Estado |
|---|---|---|
| 1 | `protein/orf.py` (6 pautas), `protein/classify.py` (4 casos), `alignment.py` para proteínas (BLOSUM62) | hecho |
| 2 | `protein/pdb.py`, `protein/hydrogens.py`, `scene/protein.py`, `protein/clients/afdb.py`, `protein/cache.py` | hecho |
| 3 | `protein/clients/blast.py` (EBI), `protein/pipeline.py`, opción **8** de la consola | hecho |
| 4 | `protein/predictors/` (ESMFold), `protein/superpose.py`, `protein/compare.py`: tu proteína predicha y comparada con la natural | hecho |
| 5a | `microscope/levels/protein.py`, motor con pantalla partida y superpuestas, `protein/session.py` | hecho |
| 5b | Web (`dnainvent_web`): herramienta Proteína 3D, trabajos asíncronos y microscopio web con proteínas | hecho |

En la consola, **8] Proteína 3D**: escribes ADN y te dice qué proteína es, cuánto se parece a
la real más cercana y la confianza de su estructura. Si no es idéntica a una real, predice la
tuya con **ESMFold** y la compara con la natural (mutaciones como `Glu7Val`, RMSD y zonas que
más cambian). Las dos estructuras se pueden guardar en PDB con hidrógenos, la tuya ya
superpuesta sobre la natural.

La búsqueda usa **BLAST del EBI**, que pide un correo electrónico real en cada búsqueda: se
usa el de Relative Developers, `relativedevelopers@gmail.com` (la variable de entorno
`DNAINVENT_EMAIL` lo sustituye). Solo se envía al EBI.

```bash
python tools/protein_info.py P68871 --save hbb.pdb     # descarga, añade H y guarda
python -m microscope --pdb hbb.pdb                     # cualquier proteína en el microscopio
```

**En el microscopio** (al final de la opción 8, o con `--pdb`): el nivel 1 es la proteína, y
desde cualquier átomo se baja igual que en el ADN (átomo → núcleo → nucleón → quarks).
En modo didáctico, la línea gruesa une los carbonos α con el color de cada tipo de
aminoácido (hidrófobo, polar, +, −, Gly/Pro); en modo realista, si es una predicción, las
nubes llevan el color de su confianza (pLDDT). Las mutaciones se marcan en magenta.
Al comparar, la natural y la tuya se ven en **pantalla partida** con la misma cámara (tecla
**V** para cambiar el lado que manejas; dentro de un átomo se ve a pantalla completa) o
**superpuestas** (tecla **C**): donde coinciden se ve blanco, donde cambian, cian y naranja.
El microscopio enseña hasta 150 aminoácidos a la vez (`scene/limits.py`).

**Hidrógenos:** AlphaFold y ESMFold no los dan; se calculan con la geometría de cada
aminoácido (`protein/templates/amino_acids.json`, generada con PyMOL). pH 7 (Lys+, Arg+, Asp-,
Glu-, extremo N -NH3+), histidina con el H en Nε y grupos que giran en posición alternada.
Comprobado con la ubiquitina 1D3Z (RMN): los H de posición fija quedan a 0,05 Å de media de los
reales. Las descargas se guardan en `~/.dnainvent/cache` (o en `DNAINVENT_CACHE`).

Datos de AlphaFold DB: CC-BY 4.0 (Jumper et al. 2021, *Nature*; Varadi et al. 2024, *NAR*).

## Base de datos de genes reales

Opción **9** de la consola: genes, variantes y enfermedades reales, guardados en local para
consultarlos y compararlos.

| Opción | Qué hace |
|---|---|
| 1 · Importar un gen | Descarga el gen de **NCBI Gene**, su ARNm de referencia y su proteína de **RefSeq** y, si es humano, sus variantes de **ClinVar** con lo que hace cada una en cada enfermedad |
| 2 · Ver un gen | Cromosoma, secuencias, variantes por clasificación y enfermedades |
| 3 · Sana frente a variante | Sitúa la variante en el ARNm (comprobando que la base de referencia es la que dice ClinVar) y enseña los codones y aminoácidos sanos y con la variante, como `GAG -> GTG`, `Glu -> Val` en la anemia falciforme |
| 4 · Buscar | Genes y enfermedades por palabras (prefijos, sin acentos); de cada enfermedad, los genes que tienen variantes que la causan |
| 5 · Genes guardados | Lista y ruta del archivo |

```
dnainvent/database/
  schema.sql      tablas: genes, secuencias, variantes, enfermedades, organismos, fármacos y dianas
  __init__.py     open_database(): crea el archivo o lo actualiza a la versión actual del esquema
  eutils.py       E-utilities del NCBI: 3 peticiones/s (10 con NCBI_API_KEY), reintentos
  ncbi.py         gen + ARNm (RefSeq Select) + proteína; secuencias bajo demanda
  clinvar.py      variantes de ClinVar y su clasificación en cada enfermedad (registros RCV)
  variants.py     situar una sustitución c.N X>Y y compararla con la secuencia sana
  queries.py      las consultas de la consola (y, más adelante, de la web)
```

**Metadatos en local, secuencias bajo demanda:** se guarda poco (los datos de cada gen y
variante) y la secuencia de un ARNm se descarga la primera vez que hace falta. El archivo está
en `~/.dnainvent/dnainvent.db` (o en `DNAINVENT_DB`) y borrarlo es seguro.

**Cada enfermedad con su clasificación:** la clasificación global de una variante en ClinVar
junta lo que han enviado todos los laboratorios. HbS es «patogénica» en global, pero para la
malaria es **protectora**. Por eso el efecto de cada variante en cada enfermedad (la causa /
factor de riesgo / protege) se lee de su registro RCV, en el XML VCV de ClinVar. Si no se puede
leer, la enfermedad se enlaza sin efecto: se sabe que se estudió, no qué hace.

**Variantes:** de momento solo se sitúan en la secuencia las sustituciones de una base dentro
de la región codificante (`c.20A>T`), que son la mayoría. Las grandes deleciones de cromosoma y
las de genes vecinos se cuentan al importar, pero no se guardan.

**Límites del NCBI:** 3 peticiones por segundo sin API key y 10 con ella. La key es gratuita
(cuenta del NCBI → *Settings → API Key Management*) y se guarda en la variable de entorno
`NCBI_API_KEY`, nunca en el código. Las peticiones llevan el mismo correo que BLAST.

Datos de NCBI Gene, RefSeq y ClinVar: dominio público. Los nombres de las enfermedades vienen
de ClinVar, en inglés.

## Usarlo desde otro proyecto

```bash
pip install -e ../DNAinvent        # en tu ordenador: los cambios del núcleo se ven al momento
```

```python
from dnainvent import sequences
sequences.transcribe("ATGC")

from dnainvent.scene.molecule import build_molecule      # necesita NumPy y AtomicWorld
build_molecule("ATGCGTACGTTAGC").n_atoms                 # 889
```

La versión del paquete es la de `dnainvent/__init__.py` (`__version__`).

## Relative Developers

DNAinvent es una app de [Relative Developers](https://relativedevelopers.com/dnainvent/).
`relativedevelopers.json` es su ficha (versión, estado, plataformas, icono y enlaces); la web
de la marca la lee con `python main.py --sync`. Cuando la web de DNAinvent tenga dirección
pública, ponla en `links.web`.

## Créditos

El modelo 3D de la doble hélice (`dnainvent/helix3d.py` y `web/static/helix3d.js` en dnainvent_web) se
inspira en la idea de [codegiovanni/DNA_helix](https://github.com/codegiovanni/DNA_helix): describir
cada hebra con una hélice paramétrica y girarla con una matriz de rotación. El código es propio (ese
repositorio no declara licencia) y usa las medidas reales del ADN B: 10,5 pares de bases por vuelta,
3,38 Å por par, 10 Å de radio y un surco mayor y otro menor.

El microscopio 3D viene del proyecto *Microscopio de ADN* (DNAmicroscope), y comparte la física
atómica y el renderizado con el visor de AtomicWorld.
