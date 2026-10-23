# Supplementary data audit: Lv et al. (2026)

**Local source:** `references/lv_2026/`  
**Article:** Lv et al. (2026), *Salinity gradients alter root-zone soil microbiome structure and nitrogen-related functional potential in alfalfa (Medicago sativa L.): a pot experiment*  
**DOI:** `10.3389/fpls.2026.1753229`

This audit covers the ten supplementary files stored with the local article copy. Dimensions below include the header row and header column where present. All `.xlsx` workbooks contain one worksheet named `Sheet1`. The accompanying `:Zone.Identifier` files are operating-system metadata and contain no scientific data.

## Supplementary file inventory

| File | Type and structure | Data level and sample identifiers | Contents and classification |
|---|---|---|---|
| `data sheet 1.csv` | CSV; 10 × 9 | Sample-level; CK-1 through MS-3 are explicit | Eight variables for each of the nine samples: EC, pH, total N, nitrate N, ammonium N, plant height, fresh weight and dry weight. Contains **soil chemistry** and **plant traits**. |
| `data sheet 2.csv` | UTF-8 CSV; 10 × 12 | Sample-level; CK-1 through MS-3 are explicit | Relative abundances for the ten displayed genera plus `Others`. Contains the **taxonomic-abundance summary** used for the dominant-genus presentation. |
| `table 1.xlsx` | Excel; `Sheet1`; 10 × 6 | Sample-level; CK-1 through MS-3 are explicit | EC (µS cm⁻¹), pH, total N (g kg⁻¹), nitrate N (mg kg⁻¹) and ammonium N (mg kg⁻¹) for all samples. This is the soil subset of `data sheet 1.csv`, with units supplied. Contains **soil chemistry**. |
| `table 2.xlsx` | Excel; `Sheet1`; 4 × 4 | Treatment-level means; rows are CK, LS and MS; individual sample IDs are absent | Mean plant height, fresh weight and dry weight for each treatment. Contains **plant-trait summaries**, not replicate-level data. |
| `table 3.xlsx` | Excel; `Sheet1`; 12 × 10 | Sample-level values, but columns are numbered 1–9 and CK-1 through MS-3 are not explicit | Transposed dominant-genus table: ten genera plus `Others` across nine unnamed columns. Contains **taxonomic-abundance summaries**. Column order is reasonably interpretable as CK-1, CK-2, CK-3, LS-1, LS-2, LS-3, MS-1, MS-2 and MS-3, but that mapping is not stated in the workbook. |
| `table 4.xlsx` | Excel; `Sheet1`; 15,806 × 10 | Sample-level matrix; columns are numbered 1–9 and explicit sample IDs are absent | 15,805 genus labels with quantitative values for nine samples. Column totals are approximately 954,000–959,000, and percentages derived from this matrix reproduce the dominant-genus values. Contains the **full genus-level taxonomic abundance matrix** available in the supplement. |
| `table 5.xlsx` | Excel; `Sheet1`; used table 28 × 3 | Replicate-level values grouped as CK, LS and MS; CK-1 through MS-3 are not explicit | Nine values each for `genus`, `kegg.entry` and `species`, organized as three values per treatment. Contains **alpha-diversity output**, but the diversity metric is not named in the file. The magnitude is compatible with the manuscript's Shannon panels, but this interpretation is not established by the workbook itself. No beta-diversity matrix or coordinates are included. |
| `table 6.xlsx` | Excel; `Sheet1`; 101 × 10 | Sample-level; CK-1 through MS-3 are explicit | Abundances for 100 entries identified only by `DBid`, including hash-like and genome/gene identifiers. By its placement relative to the manuscript's nitrogen-cycle figure, this appears to be an underlying **nitrogen-related gene/function subset**, but the file provides no functional names, process categories, KO identifiers or mapping table. It should therefore not be treated as independently annotated nitrogen-cycle data. |
| `table 7.xls` | Tab-delimited ASCII text despite the `.xls` extension; 31 × 10 | Sample-level; CK-1 through MS-3 are explicit | Thirty named metabolic pathways across nine samples, including central metabolic, biosynthetic, transport and signalling pathways. Contains **functional/KEGG pathway output**, corresponding to a selected pathway-level summary rather than a complete KO or gene catalogue. It has no workbook worksheets. |
| `table 8.docx` | Word document; 15 non-empty paragraphs; no tables | Methodological summary; no sample IDs | Defines gene-abundance terms and four equations for combining uniquely and multiply mapped reads with gene-length normalization. Contains **methodological information** only. |

No supplementary file contains a beta-diversity distance matrix, PCoA coordinates, PERMANOVA table, inferential-test table, adjusted *p*-values, assembly statistics, quality-control statistics, non-redundant gene sequences, KO-by-sample matrix or complete annotated gene catalogue.

## Environmental and plant variables

Exact replicate-level values are present in `data sheet 1.csv` for EC, pH, total N, nitrate N, ammonium N, plant height, fresh weight and dry weight. `table 1.xlsx` duplicates the five soil variables and supplies their units. `table 2.xlsx` contains treatment means for the three plant variables and is therefore a summary table rather than the primary replicate-level source.

No per-sample total-P values are present. Total P appears only as a shared baseline measurement in the main manuscript. No additional sample-level environmental or plant measurements were found beyond the eight variables in `data sheet 1.csv`.

## Methodological audit

### Functional workflow

The main manuscript reports the following sequence:

1. quality-controlled reads were assembled separately by sample with MEGAHIT, using `--min-contig-len 500`;
2. assembled sequences were described as undergoing “ORF prediction using CD-HIT” with `-n 9 -c 0.95 -G 0 -M 0 -d 0 -aS 0.9 -r 1`;
3. the longest sequence in each cluster was retained as the representative of a non-redundant gene set;
4. reads were quantified against that catalogue with Salmon;
5. representative sequences were searched against KEGG using BLASTX, after which KO results were mapped to pathways.

The supplementary methods in `table 8.docx` do **not** identify an ORF- or gene-prediction program. They clarify only the abundance calculation:

- `Ab(S) = Ab(U) + Ab(M)`;
- unique and allocated multi-mapped contributions are normalized by gene length;
- multi-mapped reads are distributed among candidate genes using a coefficient proportional to unique-mapping support.

The document does not explain how ORFs were called, how predicted genes were extracted, or whether any program preceded CD-HIT. The CD-HIT parameters and the statement that the longest sequence from each cluster was retained describe a clustering/non-redundancy step, while the manuscript labels the same step as ORF prediction. The actual ORF caller is therefore an unresolved reporting gap and must not be inferred.

The supplementary methods also do not provide a MEGAHIT version or additional assembly parameters; a Salmon version, index construction, library handling or command-line parameters; BLASTX thresholds; the KEGG database release; KO-assignment criteria; or a mapping from the `DBid` entries in `table 6.xlsx` to nitrogen-cycle processes. One line in `table 8.docx` retains the editorial instruction “state clearly: unique or multi-mapped, depending on the equation,” indicating that the definition of `M` was not fully finalized.

### Taxonomic workflow and statistics

The main manuscript reports adapter and low-quality read removal with Trimmomatic but gives no version or trimming thresholds. Host depletion used BWA-MEM v0.7.17 with `-M -k 32 -t 16`; the alfalfa reference genome accession/version is not reported. Kraken2 classification used NCBI RefSeq release 90, described as including bacterial, archaeal, fungal, viral, protist and algal genomes. Bracken was used to refine abundances from domain to species. No Kraken2 or Bracken version, database-build command, confidence threshold, read-length setting or other classification parameters are supplied.

The manuscript describes Shannon and Simpson alpha diversity at genus, species and KEGG-entry levels; Bray–Curtis PCoA and PERMANOVA for beta diversity; one-way ANOVA followed by Tukey's HSD for soil, plant and alpha-diversity comparisons; and Kruskal–Wallis tests followed by Dunn's tests with Benjamini–Hochberg correction for metagenome-derived abundance data. R 4.4.3 is reported. The number of PERMANOVA permutations, model formula and detailed implementation are absent. The supplementary files provide underlying values for only part of this workflow and do not provide formal statistical-test outputs.

## Consistency checks and unresolved issues

The following points are directly verified from the local files:

- `table 1.xlsx` exactly duplicates the five soil columns in `data sheet 1.csv`.
- CK and LS plant means in `table 2.xlsx` agree with means calculated from `data sheet 1.csv` after rounding to two decimals.
- The MS dry-weight mean does not agree: the three sample values in `data sheet 1.csv` average approximately 4.52 g, whereas `table 2.xlsx` reports 4.89 g. The source of this discrepancy is unresolved.
- `data sheet 2.csv` and percentages derived from `table 4.xlsx` agree for the displayed genera. `table 3.xlsx` also agrees except for *Sediminibacterium* in CK-3: `data sheet 2.csv` and `table 4.xlsx` support 3.409446463%, while `table 3.xlsx` contains 0.4094464625. This appears to be a transcription error, but the audit does not select a corrected value.
- The numbered columns 1–9 in `table 3.xlsx` and `table 4.xlsx` are not self-describing. Their order is strongly supported as CK-1 through MS-3 by correspondence with the explicitly labelled CSV and the derived percentages, but the workbook itself does not state that mapping.
- `table 5.xlsx` does not identify its three replicate values within each treatment or name the diversity metric.
- `table 6.xlsx` does not define `DBid` or annotate its 100 rows, preventing independent confirmation of their nitrogen-cycle roles.
- `table 7.xls` is plain tab-delimited text, not a binary Excel workbook.

## Use as full-data reference outputs

These files provide published reference values for the completed analyses or for further comparisons, subject to their stated limitations:

- **`table 4.xlsx` — principal taxonomic reference.** It provides the broadest full-depth genus-by-sample matrix in the supplement. It can support comparison of recovered genera, per-sample abundance profiles, dominant-taxon rankings and concordance with the published full-data result. The numbered sample columns must first be documented using the supported CK-1-to-MS-3 order.
- **`data sheet 2.csv` and `table 3.xlsx` — dominant-genus summaries.** They can validate the reported top-ten-plus-`Others` composition. The CSV is preferable because it carries explicit sample IDs and avoids the CK-3 *Sediminibacterium* discrepancy in `table 3.xlsx`.
- **`table 5.xlsx` — alpha-diversity benchmark with qualifications.** Its values can be compared with full-data diversity outputs only after the metric and replicate ordering are resolved. It does not provide a beta-diversity reference.
- **`table 6.xlsx` — limited gene/function reference.** It can support numerical comparison for the same 100 `DBid` entries across explicitly named samples if an authoritative identifier/function mapping is recovered. Without that mapping, it cannot validate nitrogen-process annotation.
- **`table 7.xls` — selected KEGG-pathway reference.** It can support sample-wise comparison of the 30 named pathway abundances represented in the published full-depth output. It is not a complete KEGG, KO or gene-abundance matrix.
- **`data sheet 1.csv` — covariate reference.** It is the primary sample-level source for joining soil and plant measurements to the community profiles. It is not a metagenomic reference output.

The supplement provides a full genus matrix and selected diversity, gene/function and pathway summaries, but not the complete intermediate functional outputs needed to reconstruct the published workflow exactly. The genus matrix anchors both the nine-sample ecological analysis and the CK-1/MS-1 published-versus-reconstructed comparisons; the 30-pathway table supports the functional sensitivity analysis in Report 1.
