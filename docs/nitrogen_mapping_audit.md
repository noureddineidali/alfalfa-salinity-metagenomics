# Audit of nitrogen-related DBid recoverability

## Scope and approach

Supplementary `table 6.xlsx` was read without modification. Its single worksheet contains 100 unique DBid rows and nine abundance columns (`CK-1` through `MS-3`). This audit tests only whether each DBid can be mapped automatically and reproducibly to an authoritative public record. It does not infer annotations from identifier appearance or abundance patterns.

Bare integer identifiers were tested against NCBI Gene, but matching records were retained as candidate-only because the supplement supplies no database namespace. Composite identifiers of the form `numeric prefix.locus tag` were interpreted syntactically as a candidate NCBI Taxonomy ID plus locus tag, then counted as resolved only when the locus tag and taxon matched one unique NCBI Protein record. Resolving the taxon prefix alone did not resolve the full DBid. Direct KEGG KO links were attempted for recovered NCBI protein accessions. Opaque 32-character hexadecimal strings were retained as unresolved because the supplementary material supplies no namespace or cross-reference table.

Authoritative services used: [NCBI Entrez E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/), [NCBI Gene](https://www.ncbi.nlm.nih.gov/gene/), [NCBI Protein](https://www.ncbi.nlm.nih.gov/protein/), [NCBI Taxonomy](https://www.ncbi.nlm.nih.gov/taxonomy), and the [KEGG REST API](https://www.kegg.jp/kegg/rest/keggapi.html). The per-row accessions and URLs are recorded in the mapping table.

## Recovery summary

- Resolved full DBids: **24/100 (24.0%)**
- Rows with a recovered gene name and/or protein/function description: **28/100 (28.0%)**; this comprises **24** reliable full-identifier mappings and **4** candidate-only bare-integer matches
- Rows with a direct KO mapping: **0/100 (0.0%)**
- Rows clearly nitrogen-related from the recovered annotation: **18/100 (18.0%)**
- Unresolved full DBids: **76/100 (76.0%)**

Identifier classes:

| Identifier class | Rows | Recoverability |
|---|---:|---|
| NCBI Taxonomy ID plus locus tag | 70 | 24/70 full identifiers resolved |
| bare numeric identifier; candidate NCBI Gene ID | 4 | 0/4 full identifiers resolved |
| opaque 32-character hexadecimal identifier | 26 | 0/26 full identifiers resolved |

## Comparison with Figure 7 categories

The category comparison below is intentionally strict. A category is recorded only when the recovered gene/protein annotation explicitly supports it. Generic nitrate- or nitrite-reductase wording is nitrogen-related but is not sufficient by itself to distinguish assimilatory reduction, dissimilatory reduction, or denitrification. `Organic degradation and synthesis` is assigned only to explicit glutamine-synthetase annotations; the broader contents of that paper-defined category cannot be reconstructed. `Auxiliary` is not assigned automatically because the paper does not publish a reproducible DBid-to-category rule for it.

| Figure 7 category | Directly supported DBids | Count |
|---|---|---:|
| Anammox | None | 0 |
| Nitrification | None | 0 |
| Auxiliary | None | 0 |
| Organic degradation and synthesis | 105559.Nwat_0465; 110319.CF8_1637; 1117647.M5M_14380; 1121015.N789_07300; 1122185.N792_11440 | 5 |
| Dissimilatory nitrate reduction | 110319.CF8_1167; 110319.CF8_1168 | 2 |
| Assimilatory nitrate reduction | None | 0 |
| Denitrification | 1026882.MAMP_02723; 1046627.BZARG_1503 | 2 |
| Denitrification/Dissimilatory nitrate reduction | None | 0 |
| Nitrogen fixation | 1116472.MGMO_18c00300 | 1 |

## Reproducibility assessment

The recoverable records provide evidence for individual gene or protein annotations, but they do not reconstruct the paper's complete mapping of all 100 DBids to the nine Figure 7 categories. In particular, opaque hashes cannot be linked to genes without the authors' original annotation table or gene catalogue, and broad functional categories cannot be recreated from abundance values alone. Consequently, Figure 7 cannot be reproduced faithfully from `table 6.xlsx` alone even where some individual identifiers resolve.

The mapping output preserves unresolved rows rather than substituting similarity-based guesses. A blank KO, EC number, or Figure 7 category means that no direct mapping was recovered; it does not demonstrate that the underlying sequence lacks that annotation.

## Outputs

- Mapping table: `data/processed/nitrogen_dbid_mapping.tsv`
- Reproducible resolver: `scripts/resolve_nitrogen_dbids.py`

The resolver requires the original `table 6.xlsx` at `references/lv_2026/table 6.xlsx`, internet access to the named services and Python with `openpyxl` (declared in `workflow/environment.yml`). Run it from the project root with the `alfalfa-snakemake` Python interpreter. It refreshes the mapping TSV and this audit; live database records can change, so the deposited mapping remains the reference for the published findings. Ordinary analytical reproduction uses the validated processed inputs and does not require repeating these remote queries.
