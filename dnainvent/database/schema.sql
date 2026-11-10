-- =====================================================================
-- DNAinvent — esquema de la base de datos de secuencias reales (versión 4)
-- SQLite (módulo sqlite3 de la biblioteca estándar, sin dependencias)
--
-- Idea general:
--   * En local se guardan los METADATOS (millones de filas, poco espacio).
--   * La secuencia se descarga BAJO DEMANDA y se guarda en caché
--     (columna sequences.residues, NULL hasta que alguien la abre).
--   * Las categorías salen de ontologías jerárquicas (taxonomía NCBI,
--     MONDO para enfermedades) mediante parent_id + consultas recursivas.
-- =====================================================================

PRAGMA foreign_keys = ON;

-- ---------------------------------------------------------------------
-- Fuentes y licencias: cada fila de datos sabe de dónde viene
-- ---------------------------------------------------------------------
CREATE TABLE sources (
    id              INTEGER PRIMARY KEY,
    code            TEXT NOT NULL UNIQUE,      -- 'refseq', 'ensembl', 'clinvar', 'chembl'...
    name            TEXT NOT NULL,
    url             TEXT,
    license         TEXT NOT NULL,             -- 'public-domain', 'CC-BY-4.0', 'CC-BY-SA-3.0', 'academic-only'...
    commercial_ok   INTEGER NOT NULL DEFAULT 1 CHECK (commercial_ok IN (0, 1)),
    release_version TEXT,                      -- versión de la fuente usada al importar
    imported_at     TEXT                       -- ISO 8601
);

-- ---------------------------------------------------------------------
-- Organismos (taxonomía NCBI, jerárquica)
-- ---------------------------------------------------------------------
CREATE TABLE organisms (
    taxon_id        INTEGER PRIMARY KEY,       -- ID de taxonomía NCBI (9606 = Homo sapiens)
    parent_taxon_id INTEGER REFERENCES organisms(taxon_id),
    rank            TEXT,                      -- 'species', 'genus', 'family', 'superkingdom'...
    scientific_name TEXT NOT NULL,
    name_en         TEXT,
    name_es         TEXT,
    name_ca         TEXT,
    genetic_code    INTEGER NOT NULL DEFAULT 1,   -- tabla de código genético NCBI (nuclear)
    mito_code       INTEGER                       -- código mitocondrial (2 en vertebrados)
);
CREATE INDEX idx_organisms_parent ON organisms(parent_taxon_id);

-- ---------------------------------------------------------------------
-- Genes
-- ---------------------------------------------------------------------
CREATE TABLE genes (
    id              INTEGER PRIMARY KEY,
    taxon_id        INTEGER NOT NULL REFERENCES organisms(taxon_id),
    symbol          TEXT NOT NULL,             -- 'HBB', 'TP53', 'APOE'
    name            TEXT,                      -- 'hemoglobin subunit beta'
    ncbi_gene_id    INTEGER UNIQUE,
    ensembl_id      TEXT UNIQUE,
    hgnc_id         TEXT,
    uniprot_id      TEXT,                      -- proteína revisada de UniProt: 'P68871'
    chromosome      TEXT,
    start_pos       INTEGER,                   -- coordenadas en el ensamblaje indicado
    end_pos         INTEGER,
    strand          INTEGER CHECK (strand IN (-1, 1)),
    assembly        TEXT,                      -- 'GRCh38'
    source_id       INTEGER NOT NULL REFERENCES sources(id)
);
CREATE INDEX idx_genes_symbol ON genes(symbol);
CREATE INDEX idx_genes_taxon  ON genes(taxon_id);

-- ---------------------------------------------------------------------
-- Secuencias (ADN, ARN, proteína). La secuencia en sí se baja bajo demanda.
-- ---------------------------------------------------------------------
CREATE TABLE sequences (
    id              INTEGER PRIMARY KEY,
    gene_id         INTEGER REFERENCES genes(id),          -- NULL p. ej. en genomas virales completos
    taxon_id        INTEGER NOT NULL REFERENCES organisms(taxon_id),
    accession       TEXT NOT NULL,             -- 'NM_000518', 'NP_000509', 'P68871'
    version         INTEGER,                   -- 5 en NM_000518.5 (¡las coordenadas dependen de ella!)
    mol_type        TEXT NOT NULL CHECK (mol_type IN
                        ('genomic', 'mrna', 'cds', 'ncrna', 'protein')),
    ncrna_class     TEXT,                      -- 'trna', 'rrna', 'mirna'... solo si mol_type = 'ncrna'
    is_reference    INTEGER NOT NULL DEFAULT 0 CHECK (is_reference IN (0, 1)),  -- secuencia canónica del gen
    cds_start       INTEGER,                   -- en un ARNm: dónde empieza/acaba el CDS (1-based)
    cds_end         INTEGER,
    genetic_code    INTEGER,                   -- si difiere del de organisms (mitocondrial, etc.)
    length          INTEGER,
    description     TEXT,
    residues        TEXT,                      -- NULL hasta que se descarga (caché)
    fetched_at      TEXT,
    source_id       INTEGER NOT NULL REFERENCES sources(id),
    UNIQUE (accession, version)
);
CREATE INDEX idx_sequences_gene ON sequences(gene_id);
CREATE INDEX idx_sequences_type ON sequences(mol_type);

-- Relación ARNm -> proteína (NM_000518.5 se traduce en NP_000509.1)
CREATE TABLE sequence_links (
    from_sequence_id INTEGER NOT NULL REFERENCES sequences(id),
    to_sequence_id   INTEGER NOT NULL REFERENCES sequences(id),
    relation         TEXT NOT NULL CHECK (relation IN ('translates_to', 'transcribed_from')),
    PRIMARY KEY (from_sequence_id, to_sequence_id, relation)
);

-- ---------------------------------------------------------------------
-- Variantes: SIEMPRE relativas a una secuencia concreta (accession.version)
-- ---------------------------------------------------------------------
CREATE TABLE variants (
    id              INTEGER PRIMARY KEY,
    gene_id         INTEGER NOT NULL REFERENCES genes(id),
    sequence_id     INTEGER NOT NULL REFERENCES sequences(id),   -- secuencia de referencia
    hgvs_c          TEXT,                      -- 'c.20A>T'
    hgvs_p          TEXT,                      -- 'p.Glu7Val'
    variant_type    TEXT CHECK (variant_type IN
                        ('snv', 'insertion', 'deletion', 'indel', 'duplication', 'cnv', 'other')),
    position        INTEGER,                   -- posición en la secuencia de referencia (1-based)
    ref_allele      TEXT,
    alt_allele      TEXT,
    origin          TEXT CHECK (origin IN ('germline', 'somatic', 'unknown')),  -- heredada o de tumor
    clinical_significance TEXT,                -- 'pathogenic', 'benign', 'uncertain'...
    common_name     TEXT,                      -- 'HbS', 'ΔF508', 'APOE ε4'
    dbsnp_rs        TEXT,                      -- 'rs334'
    clinvar_id      INTEGER,
    cosmic_id       TEXT,
    population_freq REAL,                      -- frecuencia alélica global (gnomAD), 0..1
    review_status   TEXT,                      -- revisión de la clasificación global (ClinVar)
    source_id       INTEGER NOT NULL REFERENCES sources(id)
);
CREATE INDEX idx_variants_gene ON variants(gene_id);
CREATE INDEX idx_variants_rs   ON variants(dbsnp_rs);

-- ---------------------------------------------------------------------
-- Enfermedades (ontología MONDO, jerárquica)
-- ---------------------------------------------------------------------
CREATE TABLE diseases (
    id              INTEGER PRIMARY KEY,
    mondo_id        TEXT UNIQUE,               -- 'MONDO:xxxxxxx'
    mesh_id         TEXT,
    parent_id       INTEGER REFERENCES diseases(id),
    name_en         TEXT NOT NULL,
    name_es         TEXT,
    name_ca         TEXT,
    category        TEXT                       -- 'cancer', 'neurodegenerative', 'blood'... (nivel alto)
);
CREATE INDEX idx_diseases_parent ON diseases(parent_id);

-- Relación variante-enfermedad. La clasificación es la de ESA enfermedad (registro RCV de
-- ClinVar), no la global de la variante: HbS es patogénica para la anemia falciforme y
-- protectora frente a la malaria. effect = NULL: se sabe que se estudió para esa
-- enfermedad, pero no hay clasificación propia (no se pudo leer su RCV).
CREATE TABLE variant_diseases (
    variant_id      INTEGER NOT NULL REFERENCES variants(id),
    disease_id      INTEGER NOT NULL REFERENCES diseases(id),
    effect          TEXT,                      -- 'causal', 'risk_factor', 'protective' o NULL
    clinical_significance TEXT,                -- 'pathogenic', 'protective'... para esta enfermedad
    review_status   TEXT,                      -- 'criteria provided, single submitter'...
    rcv_accession   TEXT,                      -- 'RCV000016574'
    source_id       INTEGER NOT NULL REFERENCES sources(id),
    PRIMARY KEY (variant_id, disease_id)
);

CREATE TABLE gene_diseases (
    gene_id         INTEGER NOT NULL REFERENCES genes(id),
    disease_id      INTEGER NOT NULL REFERENCES diseases(id),
    source_id       INTEGER NOT NULL REFERENCES sources(id),
    PRIMARY KEY (gene_id, disease_id)
);

-- ---------------------------------------------------------------------
-- Fármacos y drogas (estructura química, no secuencia; salvo péptidos)
-- ---------------------------------------------------------------------
CREATE TABLE drugs (
    id              INTEGER PRIMARY KEY,
    name            TEXT NOT NULL,
    chembl_id       TEXT UNIQUE,
    pubchem_cid     INTEGER,
    drug_type       TEXT CHECK (drug_type IN
                        ('small_molecule', 'peptide', 'protein', 'antibody',
                         'oligonucleotide', 'mrna', 'other')),
    smiles          TEXT,
    inchikey        TEXT,
    max_phase       INTEGER,                   -- 4 = aprobado (ChEMBL); 0.5 = fase 1 temprana
    sequence_id     INTEGER REFERENCES sequences(id),   -- solo fármacos peptídicos/proteicos/ARNm
    source_id       INTEGER NOT NULL REFERENCES sources(id),
    molecule_type   TEXT,                      -- tal como lo da ChEMBL: 'Small molecule', 'Cell'...
    first_approval  INTEGER,                   -- año de la primera aprobación
    withdrawn       INTEGER NOT NULL DEFAULT 0 -- 1 si se retiró del mercado (en algún país)
);
CREATE INDEX idx_drugs_name ON drugs(name);

-- Dianas: qué fármaco actúa sobre qué gen/proteína
CREATE TABLE drug_targets (
    drug_id         INTEGER NOT NULL REFERENCES drugs(id),
    gene_id         INTEGER NOT NULL REFERENCES genes(id),
    action          TEXT,                      -- 'inhibitor', 'agonist', 'antagonist', 'activator'...
    source_id       INTEGER NOT NULL REFERENCES sources(id),
    mechanism       TEXT,                      -- 'Hemoglobin S polymerisation inhibitor'
    PRIMARY KEY (drug_id, gene_id)
);

-- Mecanismos de acción tal como los da ChEMBL, aunque la diana no sea un gen guardado
-- (la hidroxiurea trata la anemia falciforme actuando sobre la ribonucleótido reductasa)
CREATE TABLE drug_mechanisms (
    drug_id         INTEGER NOT NULL REFERENCES drugs(id),
    target_chembl_id TEXT NOT NULL DEFAULT '',
    target_name     TEXT,
    target_type     TEXT,                      -- 'SINGLE PROTEIN', 'PROTEIN COMPLEX'...
    action          TEXT,
    mechanism       TEXT NOT NULL DEFAULT '',
    PRIMARY KEY (drug_id, target_chembl_id, mechanism)
);

CREATE TABLE drug_diseases (
    drug_id         INTEGER NOT NULL REFERENCES drugs(id),
    disease_id      INTEGER NOT NULL REFERENCES diseases(id),
    max_phase       INTEGER,                   -- fase máxima para esa indicación
    source_id       INTEGER NOT NULL REFERENCES sources(id),
    indication      TEXT,                      -- nombre de la indicación en ChEMBL
    efo_id          TEXT,                      -- 'MONDO:0011382', 'EFO:0000...'
    evidence        TEXT,                      -- 'label' (ficha técnica), 'trials' (solo ensayos)
    PRIMARY KEY (drug_id, disease_id)
);

-- ---------------------------------------------------------------------
-- Etiquetas libres para lo que no cubren las ontologías
-- ('antimicrobial_peptide', 'antibiotic_resistance', 'bacteriophage',
--  'vaccine_antigen', 'featured' para el catálogo didáctico...)
-- ---------------------------------------------------------------------
CREATE TABLE tags (
    id              INTEGER PRIMARY KEY,
    code            TEXT NOT NULL UNIQUE,
    name_en         TEXT,
    name_es         TEXT,
    name_ca         TEXT
);

CREATE TABLE sequence_tags (
    sequence_id     INTEGER NOT NULL REFERENCES sequences(id),
    tag_id          INTEGER NOT NULL REFERENCES tags(id),
    PRIMARY KEY (sequence_id, tag_id)
);

-- ---------------------------------------------------------------------
-- Búsqueda de texto (FTS5 viene incluido en el sqlite3 de Python)
-- ---------------------------------------------------------------------
CREATE VIRTUAL TABLE search_index USING fts5(
    kind,           -- 'gene', 'disease', 'drug', 'organism'
    ref_id UNINDEXED,
    text,
    tokenize = 'unicode61 remove_diacritics 2'
);

PRAGMA user_version = 4;
