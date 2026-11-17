"""Datos escritos a mano para los fármacos de los genes de serie, cuando lo que dan las fuentes
no basta o puede confundir.

ChEMBL marca como «retirado» (withdrawn) un fármaco que se haya retirado en CUALQUIER país, sin
decir dónde. Para los que conocemos, aquí está el alcance real, con su fuente. Para el resto,
las fichas dicen «retirado en algún país».

También corrige lo que ChEMBL no enlaza bien (fase 1d):

* APPROVALS: fármacos aprobados para una enfermedad que ChEMBL solo tiene «en ensayo» para
  ella (el elexacaftor va en una combinación, Trikafta, y ChEMBL no lo da como aprobado solo).
* CAUSE: fármacos que actúan sobre la causa aunque no actúen sobre la proteína del gen (los
  que hacen saltar un exón del ARNm de DMD, o los que hacen que SMN2 sustituya a SMN1).
* DISEASE_EQUIVALENTS: enfermedades de ClinVar que ChEMBL registra con otro nombre (ClinVar
  dice «beta Thalassemia» y ChEMBL «Thalassemia»): su ficha enseña también esos fármacos.

Los textos van en texts/*.json; aquí solo los datos.
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
    # Betibeglogén autotemcel (Zynteglo): la Comisión Europea retiró su autorización el 24-03-2022
    # a petición del fabricante, por motivos comerciales (no de seguridad). Aprobado en EE. UU.
    # en 2022. https://en.wikipedia.org/wiki/Betibeglogene_autotemcel
    "CHEMBL4594311": {"scope": "some", "year": 2022, "where": "eu"},
    # Ataluren (Translarna): autorización condicional en la UE desde 2014; la Comisión Europea
    # no la renovó y caducó el 28-03-2025 (no se confirmó que funcionara). Nunca se aprobó en
    # EE. UU. https://www.ema.europa.eu/en/medicines/human/EPAR/translarna
    "CHEMBL256997": {"scope": "some", "year": 2025, "where": "eu"},
    # Fidanacogén elaparvovec (Beqvez): Pfizer dejó de fabricarlo y venderlo en todo el mundo en
    # febrero de 2025, por motivos comerciales.
    # https://www.aabb.org/news-resources/news/article/2025/02/25/pfizer-halts-production-of-fda-approved-hemophilia-b-gene-therapy
    "CHEMBL3990034": {"scope": "worldwide", "year": 2025},
    # Valoctocogén roxaparvovec (Roctavian): BioMarin lo retiró en 2026 (EE. UU., Italia y
    # Alemania; hasta mayo de 2026) por motivos comerciales, no de seguridad.
    # https://www.bleeding.org/news/biomarin-pulls-roctavian-hemophilia-a-gene-therapy-from-us-market
    "CHEMBL4297579": {"scope": "some", "year": 2026},
}

# chembl_id -> {nombre en inglés de la enfermedad (como en la base de datos): año}
APPROVALS = {
    # Elexacaftor: aprobado en 2019 dentro de Trikafta (con tezacaftor e ivacaftor).
    # https://investors.vrtx.com/news-releases/news-release-details/fda-approves-trikafta-elexacaftortezacaftorivacaftor-and
    "CHEMBL4298128": {"Cystic fibrosis": 2019},
    # Vanzacaftor: aprobado en diciembre de 2024 dentro de Alyftrek.
    # https://investors.vrtx.com/news-releases/news-release-details/vertex-announces-us-fda-approval-alyftrektm-once-daily-next
    "CHEMBL5314934": {"Cystic fibrosis": 2024},
    # Betibeglogén autotemcel (Zynteglo): talasemia beta que necesita transfusiones (EE. UU., 2022).
    "CHEMBL4594311": {"Beta-thalassemia major": 2022},
    # Delandistrogén moxeparvovec (Elevidys): distrofia de Duchenne (EE. UU., 2023).
    "CHEMBL4650241": {"Duchenne muscular dystrophy": 2023},
    # Olaparib y talazoparib: aprobados para cánceres de personas con una mutación de BRCA1 o
    # BRCA2 heredada (ovario 2014; mama 2018).
    "CHEMBL521686": {"Hereditary breast ovarian cancer syndrome": 2014},
    "CHEMBL3137320": {"Hereditary breast ovarian cancer syndrome": 2018,
                      "Malignant tumor of breast": 2018},
}

# Actúan sobre la causa aunque su diana no sea la proteína del gen
CAUSE = {
    "CHEMBL2108278",   # eteplirsen: salto del exón 51 de DMD
    "CHEMBL4297762",   # golodirsen: exón 53
    "CHEMBL4298062",   # viltolarsen: exón 53
    "CHEMBL4297566",   # casimersen: exón 45
    "CHEMBL256997",    # ataluren: lee a través de un codón de parada prematuro
    "CHEMBL3301572",   # nusinersen: SMN2 fabrica la proteína SMN completa
    "CHEMBL4297528",   # risdiplam: ídem
    "CHEMBL3833393",   # emicizumab: hace el trabajo del factor VIII que falta
}

# nombre de ClinVar -> nombres de la base de datos con los que ChEMBL tiene sus fármacos
_BETA_THALASSEMIA = ["Thalassemia", "Beta-thalassemia major", "Beta thalassemia intermedia"]
_DUCHENNE_BECKER = ["Duchenne muscular dystrophy", "Becker muscular dystrophy"]
_SMA = ["Spinal muscular atrophy"]
DISEASE_EQUIVALENTS = {
    "beta Thalassemia": _BETA_THALASSEMIA,
    "Beta-thalassemia HBB/LCRB": _BETA_THALASSEMIA,
    "Dominant beta-thalassemia": _BETA_THALASSEMIA,
    "Beta zero thalassemia": _BETA_THALASSEMIA,
    "BETA-PLUS-THALASSEMIA": _BETA_THALASSEMIA,
    "Duchenne and Becker muscular dystrophy": _DUCHENNE_BECKER,
    "Duchenne or Becker muscular dystrophy": _DUCHENNE_BECKER,
    "Dystrophinopathies": _DUCHENNE_BECKER,
    "Becker muscular dystrophy, Cardiomyopathy, Duchenne muscular dystrophy, Dystrophin deficiency":
        _DUCHENNE_BECKER,
    "Factor VIII deficiency": ["Hereditary factor VIII deficiency disease"],
    "Severe hemophilia A": ["Hereditary factor VIII deficiency disease"],
    "Mild hemophilia A": ["Hereditary factor VIII deficiency disease"],
    "Werdnig-Hoffmann disease": _SMA,
    "Kugelberg-Welander disease": _SMA,
    "Spinal muscular atrophy, type II": _SMA,
    "Spinal muscular atrophy, type IV": _SMA,
    "Hypercholesterolemia, familial, 1": ["Familial hypercholesterolemia"],
    "LDLR-related familial hypercholesterolemia": ["Familial hypercholesterolemia"],
    "Hyperphenylalaninemia": ["Phenylketonuria"],
}

# chembl_id -> estructura del PDB con el fármaco dentro de su proteína, y el código del
# fármaco dentro de ella (para señalarlo entre las demás moléculas)
STRUCTURES = {
    "CHEMBL4101807": {"pdb": "5E83", "ligand": "5L7"},   # voxelotor (GBT440) en la hemoglobina S
}
