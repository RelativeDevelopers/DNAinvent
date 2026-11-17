"""Datos escritos a mano para los fármacos de los genes de serie, cuando lo que dan las fuentes
no basta o puede confundir.

ChEMBL marca como «retirado» (withdrawn) un fármaco que se haya retirado en CUALQUIER país, sin
decir dónde. Para los que conocemos, aquí está el alcance real, con su fuente. Para el resto,
las fichas dicen «retirado en algún país».

Los textos van en i18n; aquí solo los datos. Fase 1d: se ampliará con las explicaciones de
«por qué funciona» y las estructuras del PDB de cada fármaco.
"""

# chembl_id -> {"scope": "worldwide" | "some", "year": año, "where": código de zona (opcional)}
WITHDRAWALS = {
    # Voxelotor (Oxbryta): Pfizer lo retiró de todos los mercados en septiembre de 2024 por más
    # muertes y crisis vasooclusivas en los ensayos.
    # https://www.pfizer.com/news/press-release/press-release-detail/pfizer-voluntarily-withdraws-all-lots-sickle-cell-disease
    "CHEMBL4101807": {"scope": "worldwide", "year": 2024},
    # Crizanlizumab (Adakveo): la Comisión Europea revocó su autorización en 2023.
    # https://www.novartis.com/news/european-commission-ec-adopts-decision-endorsing-chmp-recommendation-revoke-conditional-marketing-authorization-adakveo-crizanlizumab
    "CHEMBL4297734": {"scope": "some", "year": 2023, "where": "eu"},
}

# chembl_id -> estructura del PDB con el fármaco dentro de su proteína
STRUCTURES = {
    "CHEMBL4101807": "5E83",      # voxelotor en la hemoglobina
}
