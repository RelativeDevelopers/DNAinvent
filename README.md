# 🧬 DNAinvent

**Escribe una secuencia de ADN y explórala de todas las formas posibles:** transcríbela, tradúcela a proteína, edítala con CRISPR, descubre qué proteína real se le parece… y viaja dentro de ella con un **microscopio 3D**, desde la doble hélice hasta los quarks de un protón.

DNAinvent funciona en la consola y habla **català**, **castellano** e **English**.

> Una app de [Relative Developers](https://relativedevelopers.com/dnainvent/).

---

## ✨ Qué puedes hacer

### 🔤 Trabajar con secuencias
- **Doble hélice:** a partir de una hebra, obtén su complementaria y visualízala.
- **Transcripción y traducción:** del ADN al ARN, y del ARN a la cadena de aminoácidos.
- **Codificación binaria:** convierte tu secuencia a ceros y unos.
- **Alineamiento:** compara dos secuencias y mira en qué se parecen.
- **Permutaciones:** genera todas las combinaciones posibles de una secuencia.
- **Edición CRISPR/Cas9:** corta y modifica tu secuencia como lo haría Cas9.

### 🔬 Microscopio 3D
Un viaje en primera persona dentro de tu secuencia, con cuatro niveles de zoom:

| Nivel | Escala | Qué ves |
|---|---|---|
| 1 · Hélice | nanómetros | la molécula completa, con los puentes de hidrógeno en amarillo |
| 2 · Átomo | picómetros | los orbitales del átomo y sus vecinos enlazados |
| 3 · Núcleo | femtómetros | protones y neutrones |
| 4 · Nucleón | 0,1 fm | los quarks, los gluones y el "mar" de partículas |

Puedes moverte libremente o dejar que el **piloto automático** te lleve de un nivel al siguiente (o hasta un quark del tirón). Hay cuatro modos (tecla **Tab**):

| Modo | Qué ves |
|---|---|
| **Didáctico** | para entenderlo: colores por elemento, enlaces y marcadores |
| **Realista** | más fiel a la física: densidades, sin colores inventados |
| **Atómico** | cada átomo de la hélice (oxígeno, nitrógeno, hidrógeno...) dibujado como el *Átomo en 3D* de AtomicWorld: colores por capa y todas las capas visibles. Dentro de un átomo, **O** muestra sus subniveles uno a uno y **I** sus orbitales p, con los ejes x, y, z; y puedes retroceder (**S**) hasta ver el átomo entero |
| **Plástico 3D** | electrones (e⁻), protones (+), neutrones (0), quarks up (u) y down (d) y gluones (g) como bolas de plástico brillante. Los electrones giran como en el modelo de Bohr: es un modelo de juguete, y el HUD lo explica |

### 🧪 Proteína 3D
Escribe ADN y DNAinvent te dice **qué proteína es**, cuánto se parece a la proteína real más cercana y cómo es su forma plegada (con datos de AlphaFold).
- Si tu secuencia es **inventada** o tiene **mutaciones**, predice su estructura y la compara con la natural: verás qué aminoácidos cambian (por ejemplo `Glu7Val`) y qué zonas se deforman.
- Ábrela en el microscopio y compara las dos **en pantalla partida** o **superpuestas**.
- Guarda las estructuras en formato PDB para usarlas en otros programas.

### 📚 Base de datos de genes reales
Importa genes reales del NCBI con sus variantes y las enfermedades asociadas, y consúltalos sin conexión.
- **Sana frente a variante:** compara la versión sana de un gen con la mutada, codón a codón. Por ejemplo, en la anemia falciforme: `GAG → GTG`, `Glu → Val`.
- **Buscador** de genes y enfermedades por palabras.
- Para cada variante, qué hace en cada enfermedad: si la **causa**, es **factor de riesgo** o **protege** (la variante de la anemia falciforme, por ejemplo, protege frente a la malaria).
- **Fármacos** (de ChEMBL): los que actúan sobre la proteína del gen y los indicados para cada enfermedad, en qué fase están y cómo actúan. Cada enfermedad dice si tiene tratamiento aprobado, en ensayo o nada todavía. *Aprobado quiere decir que la trata, no necesariamente que la cure. Es divulgación, no consejo médico.*
- **Un fármaco dentro de su proteína:** a partir de una estructura del PDB (por ejemplo 5E83, la hemoglobina de la anemia falciforme con voxelotor), qué aminoácidos lo rodean y cómo lo sujetan: enlace covalente, puente de hidrógeno o contacto hidrófobo, con la distancia medida en el cristal.

---

## 🚀 Instalación

### 1. Lo que necesitas
- **Python 3** instalado en tu ordenador.
- El proyecto **AtomicWorld**, que aporta la física atómica del microscopio. Si solo quieres la consola, no hace falta (ver paso 3).

### 2. Descarga los proyectos
Coloca DNAinvent y AtomicWorld uno al lado del otro, así:

```
Documents/
├── project-DNAinvent/DNAinvent
└── project-AtomicWorld/AtomicWorld
```

```bash
cd Documents/project-DNAinvent
git clone https://github.com/<tu-usuario>/DNAinvent.git
```

> También puedes descargarlo como ZIP desde el botón verde **Code → Download ZIP** de GitHub.

### 3. Instala las dependencias
Entra en la carpeta de DNAinvent y ejecuta:

```bash
pip install -r requirements.txt
```

Esto instala todo, incluido AtomicWorld desde su carpeta.

**¿Solo quieres la consola, sin microscopio?** Basta con:

```bash
pip install colorama
```

### 4. ¡A jugar!

```bash
python main.py
```

La primera vez te preguntará en qué idioma quieres usarlo y lo recordará.

---

## 🎮 Cómo se usa

### La consola
```bash
python main.py              # abre el menú
python main.py --lang ca    # en catalán (también: es, en)
```

En el menú principal tienes todas las herramientas. Algunas destacadas:
- **1] Doble hélice** — al terminar, puedes abrir la hélice en el microscopio.
- **7] Microscopio 3D**
- **8] Proteína 3D**
- **9] Base de datos de genes**
- **L] Cambiar idioma**

### El microscopio directamente
```bash
python -m microscope ATGCGTACGTTAGC                  # tu secuencia
python -m microscope ATGC --lang ca --realistic      # en catalán y en modo realista
python -m microscope ATGC --mode plastic             # empezar en otro modo: didactic, realistic, atomic o plastic
python -m microscope <secuencia> --start 120 --count 60   # un tramo de una secuencia larga
python -m microscope --pdb hbb.pdb                   # cualquier proteína guardada en PDB
```

El microscopio muestra hasta **60 pares de bases** (o 150 aminoácidos) a la vez. Si tu secuencia es más larga, eliges desde dónde empezar.

### Controles del microscopio

| Tecla | Acción |
|---|---|
| **Clic** / **Esc** | capturar / soltar el ratón (Esc otra vez: salir) |
| **W A S D** | moverse |
| **Espacio** / **Ctrl** | subir / bajar |
| **Q** / **E** | girar |
| **Shift** · **Rueda** | turbo ×5 · ajustar velocidad |
| **Enter** | 🚀 piloto automático: bajar un nivel hacia lo que apuntas |
| **Mayús + Enter** | 🚀 viaje completo hasta un quark |
| **Retroceso** | subir un nivel |
| **Tab** · **Mayús + Tab** | modo: didáctico → realista → atómico → plástico 3D (y al revés) |
| **O** · **I** | (modo atómico, dentro de un átomo) un subnivel cada vez · un orbital p cada vez |
| **L** | cambiar idioma |
| **H** · **F11** | ayuda · pantalla completa |
| **C** | (comparando proteínas) pantalla partida ↔ superpuestas |
| **V** | (pantalla partida) cambiar el lado que controlas |

---

## 💡 Consejos

- **¿Tienes gráfica NVIDIA en Windows?** Para que el microscopio vaya fluido: *Configuración → Sistema → Pantalla → Gráficos*, añade el `python.exe` que usas y elige **Alto rendimiento**. Al arrancar, la consola te dice qué tarjeta gráfica está usando.
- **Base de datos más rápida:** el NCBI permite 3 consultas por segundo, o 10 si tienes una API key gratuita (cuenta del NCBI → *Settings → API Key Management*). Guárdala en la variable de entorno `NCBI_API_KEY`.
- **Tus datos se guardan en** `~/.dnainvent/`: el idioma elegido, la base de datos de genes y las descargas. Puedes borrar la base de datos o la caché cuando quieras; se vuelven a crear.

---

## 🔒 Privacidad y conexión a Internet

Casi todo funciona sin conexión. Solo se conecta a Internet para:
- **Identificar proteínas** (búsqueda BLAST en el EBI) y **descargar estructuras** de AlphaFold DB.
- **Importar genes** del NCBI y **fármacos** de ChEMBL (con UniProt para saber cuál es la proteína de cada gen).
- **Descargar estructuras con fármacos** del PDB (RCSB).

El EBI y el NCBI piden un correo de contacto en cada búsqueda; se usa el de Relative Developers (`relativedevelopers@gmail.com`), y puedes cambiarlo con la variable de entorno `DNAINVENT_EMAIL`. No se envía nada más.

---

## 🧠 Cuánto hay de real

DNAinvent intenta enseñar ciencia de verdad:
- La doble hélice se construye átomo a átomo con medidas reales del ADN forma B (**10,5 pares por vuelta**, 3,38 Å por par). Los puentes de hidrógeno no están dibujados a mano: salen de la geometría (2 en A–T, 3 en G–C).
- Los orbitales, los núcleos y los quarks usan modelos físicos reales, compartidos con **AtomicWorld**.

Y también tiene sus simplificaciones: los orbitales son los de un átomo aislado, la hélice es ideal (el ADN real varía según la secuencia), el núcleo se modela como una "gota líquida" y el movimiento de los quarks es una ilustración, porque en física cuántica no tienen trayectorias.

---

## 👩‍💻 Para desarrolladores

```bash
python -m pytest            # ejecutar los tests
pip install -e .            # instalar DNAinvent como paquete en otro proyecto
```

```python
from dnainvent import sequences
sequences.transcribe("ATGC")
```

**Estructura del proyecto:**

```
dnainvent/     lógica biológica (secuencias, proteínas, base de datos, física de la escena)
microscope/    microscopio 3D de escritorio
i18n/          traducciones (un JSON por idioma)
ui/            consola y ajustes
tests/         tests automáticos
tools/         scripts auxiliares
```

**¿Quieres añadir un idioma?** Copia `i18n/locales/es.json` con el código del nuevo idioma y traduce los textos. Los tests comprueban que no falte ninguno.

La versión web de DNAinvent está en su propio repositorio, `dnainvent_web`, que usa este como paquete.

---

## 🙏 Créditos y datos

- **Datos de proteínas:** AlphaFold DB, licencia CC-BY 4.0 (Jumper et al. 2021, *Nature*; Varadi et al. 2024, *NAR*).
- **Datos de genes:** NCBI Gene, RefSeq y ClinVar (dominio público). Los nombres de las enfermedades vienen de ClinVar, en inglés.
- **Fármacos:** ChEMBL, licencia CC BY-SA 3.0 (Zdrazil et al. 2024, *NAR*). **Proteínas de cada gen:** UniProt, CC BY 4.0. **Estructuras con fármacos:** RCSB PDB (dominio público).
- **Modelo de la doble hélice:** inspirado en la idea de [codegiovanni/DNA_helix](https://github.com/codegiovanni/DNA_helix) (hélice paramétrica girada con una matriz de rotación). El código es propio y usa las medidas reales del ADN B.
- **Microscopio 3D:** viene del proyecto *Microscopio de ADN* (DNAmicroscope) y comparte física y renderizado con el visor de **AtomicWorld**.

---

<p align="center">Hecho con 🧬 y curiosidad por <a href="https://relativedevelopers.com/dnainvent/">Relative Developers</a></p>