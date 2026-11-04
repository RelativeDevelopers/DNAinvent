"""Física de lo que dibuja el microscopio 3D (escritorio y web), sin OpenGL.

Necesita NumPy y la física atómica de AtomicWorld:

    pip install dnainvent[scene]

El resto del núcleo de DNAinvent no tiene dependencias.

    molecule.py   secuencia -> todos los átomos de la doble hélice (forma B, 10,5 pares/vuelta)
    elements.py   H, C, N, O y P: datos de la tabla periódica y su isótopo más abundante
    describe.py   resumen de la molécula para el visor web (/api/structure de dnainvent_web)

Del átomo hacia abajo (orbitales, núcleo, nucleones) la física es la de AtomicWorld
(atomicworld.scene): así las dos apps enseñan exactamente lo mismo.
"""
