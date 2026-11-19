"""El visor 3D de escritorio de «Cuerpo humano» (python -m bodyviewer): el cuerpo de referencia
pieza a pieza, con capas por sistema; al pulsar una pieza dice qué es, qué células tiene y qué
genes las marcan. Lo abre la opción A de la consola.

    scene.py      las mallas (meshes/*.dnab) listas para la tarjeta gráfica
    camera.py     la cámara que gira alrededor del cuerpo
    renderer.py   el dibujo (moderngl) y saber qué pieza hay bajo el ratón
    hud.py        los textos encima (Pillow)
    app.py        la ventana (pyglet), el ratón y las teclas
"""
