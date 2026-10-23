# Raw-read taxonomic comparison

This directory contains the reproducible CK-1 comparison between the genus profile released by Lv et al. (2026) and the independent PlusPF-8 Bracken profile. The published CK-1 column in `table 4.xlsx` was identified using the explicitly labelled `data sheet 2.csv`. Report 2 extends the descriptive comparison to MS-1; its two-sample figures are built by [`scripts/build_two_sample_raw_read_figures.py`](../../scripts/build_two_sample_raw_read_figures.py) from the validated profiles and metadata.

## Run

From the directory containing this README, with Python packages `numpy`, `scipy` and `matplotlib` available. The validation run used Python 3.13.0, NumPy 2.3.2, SciPy 1.16.1 and Matplotlib 3.10.5.

```bash
python compare_ck1_profiles.py \
  --published ck1_published_table4_extracted.csv \
  --bracken CK1_pluspf8_c0_genus.bracken \
  --outdir .
```

This rebuilds the comparison CSV tables, JSON metrics, figures and Markdown report. Inputs have been included for exact reproducibility.

## Files

- `CK1_comparison_report.md`: complete results-style narrative, methods and figure references.
- `all_genus_exact_label_comparison.csv`: union of published and independent label sets; empty counterpart cells mean label absent in the other table, **not biological absence**.
- `shared_positive_genera.csv`: matched positive-abundance labels used in correlation.
- `top20_each_profile.csv`: top-20 entries from each profile, with counterpart abundance/rank where present.
- `comparison_metrics.json`: all principal descriptive statistics.
- `fig1_...` through `fig5_...`: publication-style figures.
- `ck1_published_table4_extracted.csv`: source numeric column, values unchanged; relative percentages are calculated in script.
- `CK1_pluspf8_c0_genus.bracken`: genus-level Bracken table, retained here to make this comparison self-contained. It is byte-identical to the archived output under `results/taxonomy/raw-read-pilot/`.

## Methodology

Comparison is **one sample versus the same sample** through two distinct bioinformatic pipelines. Genus labels are matched *exactly*; taxonomic synonyms and provisional names are not collapsed. Both abundance columns are normalized to sum one within their own source. Thus Spearman and exact-label Bray–Curtis are **descriptive measures of pipeline concordance**, not treatment-level tests or evidence for biochemical function. The independent PlusPF-8 run classified 5.42% of the 33,808,912 input non-host fragments; interpreting genus proportions as whole-community cell fractions is inappropriate.

Twelve genus-label cells in the source workbook were stored as Excel date-formatted numbers. Their underlying serial values are retained in the extracted input because the intended original labels are not recoverable. A sensitivity rerun using the processed table's ISO date representations produced identical metrics because none of these 12 labels matched a Bracken genus.

Database release, MD5 validation, read-pair counts and executed commands are documented in the [raw-read provenance](../../docs/raw_read_provenance.tsv) and [Report 2](../../report/alfalfa_ck1_raw_read_reconstruction.html). The original publication does not specify its exact Kraken database build or abundance denominator, limiting exact computational replication.
