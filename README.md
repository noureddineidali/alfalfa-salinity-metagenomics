# Salinity-associated restructuring of the alfalfa root-zone microbiome

How does salinity relate to the alfalfa root-zone microbiome—and how much of a published metagenomic profile can be reconstructed from public data? This project examines both questions through complementary reanalyses of Lv et al. (2026): all nine released sample-level profiles and independent raw-read reconstructions of CK-1 and MS-1.

**[Explore Report 1 · nine-sample ecological reanalysis →](https://noureddineidali.github.io/alfalfa-salinity-metagenomics/report/alfalfa_salinity_metagenome_reanalysis.html)**

**[Explore Report 2 · CK-1 and MS-1 raw-read reconstruction →](https://noureddineidali.github.io/alfalfa-salinity-metagenomics/report/alfalfa_ck1_raw_read_reconstruction.html)**

[BioProject PRJNA1399466](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1399466) · [Publication DOI](https://doi.org/10.3389/fpls.2026.1753229)

<a href="results/figures/project_graphical_abstract.svg"><img src="results/figures/project_graphical_abstract.svg" alt="Graphical abstract showing the nine-sample published-data reanalysis alongside CK-1 and MS-1 raw-read reconstruction, their workflows, principal findings and interpretation limits" width="100%"></a>

## Two reports, one public dataset

| Report 1 · published data | Report 2 · deposited reads |
|---|---|
| Nine published sample profiles (CK, LS and MS; three per treatment): community composition, diversity and functional pathways. | CK-1 and MS-1 reconstructed from SRA reads through QC, host depletion and Kraken2/Bracken profiling; compared descriptively with published profiles. |

The sampled material was **root-zone soil**, not a strictly isolated rhizosphere fraction.

## Findings at a glance

| Genus-level community structure | Relative functional structure | Non-host fragments classified |
|:---:|:---:|:---:|
| **R² = 0.775** | **R² = 0.822** | **5.17–5.42%** |
| Nine-sample Bray–Curtis PERMANOVA | Nine-sample, composition-normalized 30-pathway PERMANOVA | CK-1 and MS-1, Kraken2 PlusPF-8 |

In the functional analysis, raw PCoA1 tracked summed selected-pathway abundance (**r = −0.9988**); the normalized result above describes relative pathway structure.

## Reproducibility

The published-data analyses use **R 4.3.3** and the synchronized [`renv.lock`](renv.lock):

```r
renv::restore()
```

The raw-read workflow keeps bioinformatics tools in [`environment.yml`](environment.yml) (`alfalfa-meta`) and Snakemake plus validation libraries in [`workflow/environment.yml`](workflow/environment.yml) (`alfalfa-snakemake`). With `alfalfa-meta` active, select a sample from the [canonical manifest](data/metadata/sample_manifest.tsv):

```bash
"$(conda info --base)/envs/alfalfa-snakemake/bin/snakemake" \
  --snakefile workflow/Snakefile \
  --profile workflow/profiles/local \
  --config 'active_samples=["MS-1"]'
```

If Conda environments are installed outside the base `envs/` directory, substitute the path to `alfalfa-snakemake/bin/snakemake`. The [workflow guide](workflow/README.md) covers setup and sample selection. For optional live monitoring, run `snakesee watch .` in a separate terminal with `alfalfa-snakemake` active.

## Explore the repository

[`analysis/`](analysis/) and [`data/processed/`](data/processed/) support Report 1; [`workflow/`](workflow/), [`scripts/`](scripts/) and [`results/`](results/) contain the raw-read workflow and outputs. The [sample manifest](data/metadata/sample_manifest.tsv) links samples to public accessions.

Inference is limited by three sequenced samples per treatment in Report 1 and one reconstructed sample per condition in Report 2. Published abundance units and nitrogen-function mappings remain incomplete; most non-host fragments were unclassified with the workstation-scale Kraken2 database. Supporting audits: [dataset](docs/dataset_reconnaissance.md), [supplementary files](docs/supplementary_data_audit.md), [processed inputs](docs/processed_data_validation.md) and [nitrogen identifiers](docs/nitrogen_mapping_audit.md).

