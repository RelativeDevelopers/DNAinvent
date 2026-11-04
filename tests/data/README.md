# Datos de prueba

| Archivo | Qué es | Origen y licencia |
|---|---|---|
| `AF-P68871-F1-model_v6.pdb` | Hemoglobina beta humana (147 aminoácidos) predicha por AlphaFold, sin hidrógenos | [AlphaFold DB](https://alphafold.ebi.ac.uk/entry/P68871), CC-BY 4.0 (Jumper et al. 2021; Varadi et al. 2024) |
| `1d3z_model1.pdb` | Ubiquitina humana medida por RMN, **con hidrógenos** (solo el modelo 1 de 10) | [PDB 1D3Z](https://www.rcsb.org/structure/1D3Z), CC0 (Cornilescu et al. 1998) |
| `clinvar_esummary_hbb.json` | Respuesta real de ClinVar (esummary JSON) para 3 variantes de HBB, de noviembre de 2026 | [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/), dominio público |
| `clinvar_vcv_msh2.xml` | Respuesta real de ClinVar (XML VCV) de la variante MSH2 c.2047G>A, con sus 4 registros RCV, de noviembre de 2026 | [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/variation/90868/), dominio público |

La ubiquitina sirve para comprobar que los hidrógenos calculados (`dnainvent/protein/hydrogens.py`)
quedan donde están los reales.
