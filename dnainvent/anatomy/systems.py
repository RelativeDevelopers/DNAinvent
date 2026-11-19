"""Los sistemas del cuerpo: las capas que se encienden y apagan en el visor.

BodyParts3D no dice a qué sistema va cada pieza; se deduce de su nombre y de sus dos árboles
(«es un» y «forma parte de»), con reglas en este orden:
    1. palabras del nombre muy claras (nervio, arteria, ligamento, «eyeball»...);
    2. la clase a la que pertenece en el árbol «es un» (músculo, hueso, arteria, vena...);
    3. lo mismo para aquello de lo que forma parte («cabeza larga del bíceps» → bíceps → músculo).
Lo que no encaja queda en «other». Los nombres de los sistemas en cada idioma están en i18n
(anatomy.systems.<clave>), no aquí.
"""

import re

# (clave, color en el modo didáctico). El orden es el de la lista de capas.
SYSTEMS = (
    ("skeletal", "#e8dcc4"),
    ("muscular", "#b8433c"),
    ("connective", "#d9d4c7"),
    ("heart", "#8e1b24"),
    ("arteries", "#d22f2f"),
    ("veins", "#3557b7"),
    ("nervous", "#e2c044"),
    ("respiratory", "#e79aa6"),
    ("digestive", "#d98a4b"),
    ("urinary", "#c9a227"),
    ("reproductive", "#c77d9b"),
    ("endocrine", "#8e6cc2"),
    ("lymphatic", "#5aa36a"),
    ("senses", "#4fa3b8"),
    ("skin", "#e9c2a6"),
    ("other", "#9a9a9a"),
)
SYSTEM_KEYS = tuple(k for k, _ in SYSTEMS)

# 1. Palabras del nombre (en minúsculas), por orden. Los vasos van primero: «arteria
#    cerebelosa» es una arteria, no cerebelo.
_NAME_RULES = (
    ("arteries", r"arter|aorta|palmar arch"),
    ("veins", r"\bvein|venous|vena\b"),
    ("heart", r"^heart$|wall of heart|mitral|tricuspid|pulmonary valve|aortic valve|(left|right) ventricle|atrium|"
              r"papillary muscle|wall of ventricle"),
    ("skeletal", r"\bbone\b|vertebra\b|\brib\b|sternum|manubrium|xiphoid|cartilage|meniscus|"
                 r"intervertebral|sacrum|coccyx|skull|mandible|maxilla|clavicle|\bscapula\b|patella"),
    ("senses", r"eyeball|choroid$|cornea|\biris\b|\bretina\b|\blens\b|sclera|ciliary body|corona ciliaris|"
               r"lacrimal (lake|sac|canaliculus|gland)|nasolacrimal|cochlea|tympan|vitreous|"
               r"conjunctiva|tarsal plate|check ligament|common tendinous ring|trochlea of"),
    ("nervous", r"nerve|ganglion|thalamus|lateral ventricle|third ventricle|fourth ventricle|"
                r"interventricular foramen|interpeduncular|tentorium|choroid plexus|cerebral hemisphere|"
                r"cerebell|brainstem|\bpons\b|medulla oblongata|spinal cord|cerebral aqueduct|"
                r"hypothalamus|hippocamp|amygdala|caudate nucleus|putamen|cerebral cortex|midbrain|colliculus"),
    ("skin", r"^skin$|\bhair\b|eyebrow|^lip$|external ear|auricle"),
    ("muscular", r"muscle|\bhead of\b|trapezius|deltoid|pectoralis|tensor|constrictor|"
                 r"pharyngeus|glossus\b|levator|sphincter"),
    ("reproductive", r"penis|testis|epididym|deferent|seminal|prostat|scrot"),
    ("respiratory", r"epiglott|bronch|\blung|trache|laryn|pleura|nasal cavity"),
    ("urinary", r"kidney|ureter|urinary bladder|urethra"),
    ("endocrine", r"pituitary|adrenal|thyroid gland|parathyroid"),
    ("lymphatic", r"spleen|thymus|lymph|tonsil"),
    ("digestive", r"tooth|gingiva|tongue|stomach|intestin|colon|rectum|esophag|pancrea|gallbladder|"
                  r"duoden|jejun|ileum|cecum|appendix|pharyn|palate|liver|hepat|sublingual|salivary|"
                  r"parotid|submandibular gland|mesentery|biliary"),
    ("connective", r"ligament|tendon|raphe|membrane|iliotibial tract|retinaculum|linea alba|fascia|"
                   r"aponeurosis|conus elasticus"),
)
# 2. Clases del árbol «es un».
_CLASS_RULES = (
    ("muscular", {"muscle organ", "head of muscle organ", "zone of muscle organ", "set of muscles", "musculature"}),
    ("skeletal", {"bone organ", "cartilage organ", "intervertebral disk"}),
    ("heart", {"heart", "cardiac chamber", "cardiac valve", "wall of heart"}),
    ("arteries", {"artery", "segment of arterial tree organ"}),
    ("veins", {"vein", "segment of venous tree organ"}),
    ("nervous", {"nerve", "segment of neural tree organ", "segment of neuraxis", "cell part cluster of neuraxis", "ganglion"}),
    ("respiratory", {"segment of tracheobronchial tree", "lung"}),
    ("digestive", {"tooth", "zone of small intestine"}),
)
_NAME_RULES = tuple((key, re.compile(rx)) for key, rx in _NAME_RULES)


def classify(bp, fma: str, _depth: int = 0) -> str:
    """El sistema de un concepto de BodyParts3D (bp: un BodyParts3D o algo con name(),
    ancestors(), isa_parents() y partof_parents())."""
    name = bp.name(fma).lower()
    for key, rx in _NAME_RULES:
        if rx.search(name):
            return key
    classes = {bp.name(a).lower() for a in bp.ancestors(fma, "isa")}
    for key, names in _CLASS_RULES:
        if classes & names:
            return key
    if _depth < 4:
        for parent in sorted(bp.partof_parents(fma)) + sorted(bp.isa_parents(fma)):
            found = classify(bp, parent, _depth + 1)
            if found != "other":
                return found
    return "other"
