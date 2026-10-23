# Processed data validation

The canonical inputs are the validated, source-derived TSV matrices in `data/processed/` and the accession manifest in `data/metadata/`. They preserve the source numerical values while making sample identifiers, column names and metadata explicit. They are processed inputs, rather than untouched copies of the supplementary files. Reproducing the published analyses starts from these TSV files; the initial workbook-conversion utility is not required.

| Input | Source and preparation |
|---|---|
| `data/metadata/sample_manifest.tsv` | BioSample/SRA accessions documented in [dataset reconnaissance](dataset_reconnaissance.md); explicit treatments and replicate numbers. |
| `data/processed/environmental_plant_metadata.tsv` | `data sheet 1.csv`; descriptive variable headers with units, plus treatment and replicate fields. |
| `data/processed/genus_abundance.tsv` | `table 4.xlsx`; numbered columns mapped using `data sheet 2.csv`, as validated below. |
| `data/processed/pathway_abundance.tsv` | Tab-delimited `table 7.xls`; first header renamed `pathway`. |
| `data/processed/nitrogen_dbid_mapping.tsv` | Identifier-recovery audit of `table 6.xlsx`, documented in the [nitrogen mapping audit](nitrogen_mapping_audit.md). |

Original files are available with [Lv et al. (2026)](https://doi.org/10.3389/fpls.2026.1753229); their formats, contents and inconsistencies are documented in the [supplementary-data audit](supplementary_data_audit.md).

## Sample structure

- `sample_manifest.tsv` contains exactly 9 unique samples: CK-1–CK-3, LS-1–LS-3, and MS-1–MS-3.
- Treatment assignments are Control (CK), Low salinity (LS), and Moderate salinity (MS), with replicates 1–3 in each treatment.
- `environmental_plant_metadata.tsv` contains the same 9 samples and 8 original measured variables. It has 9 data rows and 11 columns after adding sample, treatment, and replicate fields.

## Genus matrix

- Dimensions: 15,805 genus rows × 9 sample columns (10 columns including the genus label).
- The numbered workbook columns were verified against the ten explicitly labelled dominant genera in `data sheet 2.csv`: 1=CK-1, 2=CK-2, 3=CK-3, 4=LS-1, 5=LS-2, 6=LS-3, 7=MS-1, 8=MS-2, 9=MS-3.
- All 90 genus–sample percentage comparisons agreed within 5E-8 percentage points. The largest absolute difference was 0.000000000498 percentage points.
- Missing genus labels: 0.
- Duplicated genus labels: 0.
- Missing abundance cells: 0.

### Preserved abundance-column totals

| Sample | Column total |
|---|---:|
| CK-1 | 954879.479053 |
| CK-2 | 955206.870793 |
| CK-3 | 957597.132289 |
| LS-1 | 958532.457664 |
| LS-2 | 956293.566292 |
| LS-3 | 958147.416404 |
| MS-1 | 953627.079059 |
| MS-2 | 954192.700428 |
| MS-3 | 954064.468387 |

During input preparation, no abundance values were normalized, transformed, filtered, aggregated, or statistically analysed. Percentages were calculated only in memory to verify the sample-column mapping and were not written to the genus matrix.

## Pathway matrix

- Source format: ASCII tab-delimited text (not a binary Excel workbook).
- Dimensions: 30 pathway rows × 9 sample columns (10 columns including the pathway label).
- Column names: `pathway`, `CK-1`, `CK-2`, `CK-3`, `LS-1`, `LS-2`, `LS-3`, `MS-1`, `MS-2`, `MS-3`.
- The 30 pathway names and their original row order were preserved exactly; only the blank first-column header was renamed `pathway`.
- All 270 numerical cells match the source text exactly.
- Missing pathway names: 0.
- Duplicated pathway names: 0.
- Missing abundance cells: 0.
- Negative abundance values: 0.

Validation passed. Values were not transformed, converted to percentages, normalized, scaled, rounded, filtered, reordered, renamed, or imputed. The source does not define the abundance unit, so the values remain labelled as pathway abundance without a more specific unit.

## Regenerating analytical outputs

After `renv::restore()`, run the numbered R scripts from the project root. Scripts 02–06 export PNG figures and TSV results; 05 reads the PCoA tables exported by 04. To compare regenerated outputs without replacing the published files:

```bash
ALFALFA_RESULTS_DIR=/tmp/alfalfa-results Rscript -e \
  'for (p in sort(list.files("analysis", pattern="^0[2-6]_.*[.]R$", full.names=TRUE))) source(p)'
```

Without the override, exports go to `results/`. The final report's six-genus panel and refined taxonomic PCoA are produced by `scripts/build_report1_presentation_figures.py`; the latter uses the saved PCoA coordinates. The raw-read report's three comparison figures are produced by `scripts/build_two_sample_raw_read_figures.py`. Both Python scripts use the scientific Python packages in `workflow/environment.yml`.

Export verification reproduced all eight retained R PNGs pixel-for-pixel and all eight existing R-derived TSV tables exactly. The canonical input matrices and published result files were not overwritten during this check.
