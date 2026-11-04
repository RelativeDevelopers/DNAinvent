"""De la secuencia de ADN a la proteína plegada en 3D.

    orf        buscar marcos de lectura abiertos en las 6 pautas y traducirlos
    classify   decidir, con los números de una búsqueda BLAST, si la proteína es conocida,
               una variante, un pariente lejano o inventada
    pdb        leer y escribir archivos PDB
    hydrogens  añadir los hidrógenos que faltan (tabla templates/amino_acids.json)
    geometry   vectores y coordenadas internas, en Python puro
    confidence franjas de confianza (pLDDT) de una predicción
    cache      caché en disco (~/.dnainvent/cache)
    clients/   servicios externos: afdb (AlphaFold DB)

La molécula 3D para el microscopio (con NumPy) está en dnainvent/scene/protein.py.

Como el resto del paquete, solo contiene lógica: no pregunta ni imprime nada.
"""
