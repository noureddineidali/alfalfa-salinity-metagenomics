# Dataset reconnaissance: PRJNA1399466

**Metadata accessed:** 7 October 2026

## 1. Study overview

[PRJNA1399466](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1399466) is an NCBI BioProject titled *Shotgun metagenomics of alfalfa root-zone soil under a salinity gradient*. It is registered as a multispecies project containing raw sequence reads from a controlled pot experiment with alfalfa (*Medicago sativa* L.). The study was designed to characterize microbial community composition and functional potential across a salinity gradient, with particular attention to nitrogen-cycling functions and KEGG pathways, and to relate the metagenomic profiles to soil chemistry and plant growth measurements.

The associated publication is Lv et al. (2026), “Salinity gradients alter root-zone soil microbiome structure and nitrogen-related functional potential in alfalfa (*Medicago sativa* L.): a pot experiment,” published in *Frontiers in Plant Science* (DOI: [10.3389/fpls.2026.1753229](https://doi.org/10.3389/fpls.2026.1753229)).

The sampled material requires careful terminology. The authors define it as root-zone soil collected by coring within approximately 3 cm of visible roots. Soil tightly adhering to the root surface was not physically separated. The samples should therefore be described as operationally defined root-zone soil, not as a strictly isolated rhizosphere fraction.

## 2. Experimental design

The experiment used alfalfa cultivar ‘Algonquin’. Field-collected chestnut soil was sieved, homogenized and air-dried before being placed into pots. Each pot contained 7.5 kg of soil and received 20 seedlings. The pots were irrigated daily with 500 mL of deionized water.

Three salinity treatments were applied. CK was the non-saline control. LS received 3.0 g salt kg⁻¹ dry soil, described in the article as 0.3% salinity, and MS received 4.5 g salt kg⁻¹ dry soil, described as 0.45% salinity. The added salt was a 1:1 mixture by weight of NaCl and Na₂SO₄, mixed into the soil once at the beginning of the treatment. Treatments began on 17 June 2024 and continued until harvest on 30 September 2024, a period of 105 days.

The sequencing dataset contains three biological replicates per treatment: CK-1 to CK-3, LS-1 to LS-3 and MS-1 to MS-3. This gives nine biological samples in total. The experimental schematic in the publication depicts eight pots per treatment, whereas only three samples per treatment were sequenced and used in the reported metagenomic comparisons. The available metadata do not explain how the three sequenced samples were selected from the eight pots, nor do they explicitly map each sample to a pot identifier.

## 3. Sequencing dataset

Each BioSample is associated one-to-one with an SRA sample, experiment and run. “Read pairs” in the table correspond to SRA spots for which both mates are present. Each pair comprises two individual reads.

| Sample | BioSample | SRA sample | Experiment | Run | Treatment | Replicate | Read pairs | Bases | NCBI archive size (GB) | ENA FASTQ.gz size (GB) |
|---|---|---|---|---|---|---:|---:|---:|---:|---:|
| CK-1 | SAMN54495940 | SRS27706072 | SRX31738113 | SRR36739085 | Control | 1 | 34,283,431 | 10,285,029,300 | 5.850 | 7.640 |
| CK-2 | SAMN54495941 | SRS27706074 | SRX31738114 | SRR36739084 | Control | 2 | 36,814,679 | 11,044,403,700 | 6.335 | 8.256 |
| CK-3 | SAMN54495942 | SRS27706073 | SRX31738115 | SRR36739083 | Control | 3 | 33,785,543 | 10,135,662,900 | 5.819 | 7.592 |
| LS-1 | SAMN54495943 | SRS27706075 | SRX31738116 | SRR36739082 | Low salinity | 1 | 42,849,418 | 12,854,825,400 | 7.152 | 9.346 |
| LS-2 | SAMN54495944 | SRS27706076 | SRX31738117 | SRR36739081 | Low salinity | 2 | 45,851,914 | 13,755,574,200 | 7.821 | 10.213 |
| LS-3 | SAMN54495945 | SRS27706077 | SRX31738118 | SRR36739080 | Low salinity | 3 | 36,774,905 | 11,032,471,500 | 6.188 | 8.078 |
| MS-1 | SAMN54495946 | SRS27706078 | SRX31738119 | SRR36739079 | Moderate salinity | 1 | 34,492,089 | 10,347,626,700 | 5.779 | Not available in ENA report |
| MS-2 | SAMN54495947 | SRS27706080 | SRX31738120 | SRR36739078 | Moderate salinity | 2 | 35,455,278 | 10,636,583,400 | 5.924 | 7.758 |
| MS-3 | SAMN54495948 | SRS27706079 | SRX31738121 | SRR36739077 | Moderate salinity | 3 | 35,625,707 | 10,687,712,100 | 5.931 | 7.766 |

The complete NCBI dataset comprises nine biological samples and nine sequencing runs. It contains **335,932,964 SRA spots/read pairs**, equivalent to approximately **671,865,928 individual reads**, and a verified total sequence yield of **100,779,889,200 bases (100.78 Gb)**. The SRA study accession is `SRP660808`, and the submission accession is `SRA2312055`.

Sequencing was performed on an Illumina NovaSeq 6000 in paired-end 150 bp mode. The 300 bp average SRA spot length is consistent with two 150 bp reads per spot. The SRA records describe the libraries as WGS strategy, METAGENOMIC source, RANDOM selection and PAIRED layout. Approximately 1 µg of DNA per sample was fragmented to about 450 bp before library construction. The SRA RunInfo insert-size field is nevertheless recorded as zero and is not informative.

## 4. Sample and environmental metadata

All BioSamples are classified as `soil metagenome` (NCBI taxonomy ID 410658). The host is recorded as alfalfa, the isolation source as root-zone soil, and the collection date as 30 September 2024. The geographic location is Saihan District, Hohhot, Inner Mongolia Autonomous Region, China, at 40.807222 N, 111.746667 E. The BioSample environmental descriptors are agricultural biome (`ENVO:00000077`) at broad scale, alfalfa root-zone soil at local scale and soil (`ENVO:00001998`) as the environmental medium. The samples use the MIxS soil metagenome package `MIMS.me.soil.6.0`. Sampling depth and elevation are recorded as “not collected.”

Additional experimental metadata are available from the publication rather than from the BioSample records. These include the ‘Algonquin’ cultivar, chestnut soil classification, pot dimensions, quantity of soil and number of seedlings per pot, the salt composition and dose, treatment timing, irrigation regime and operational definition of the sampled soil compartment.

Baseline measurements shared by the pots included soil electrical conductivity, pH, nitrate nitrogen, total nitrogen and total phosphorus. At harvest, the reported soil variables were electrical conductivity, pH, nitrate nitrogen, ammonium nitrogen and total nitrogen. Plant measurements included height, fresh biomass and dry biomass. The BioSample records do not contain the per-sample values for these variables; the associated publication supplements appear to be their public source.

## 5. Public data availability

The raw sequence data are public through NCBI SRA under PRJNA1399466. NCBI provides a complete record for all nine runs. The sum of the run-level NCBI archive sizes is approximately 56.8 GB, while the BioProject summary displays approximately 59,561 MB. These figures describe archive representations rather than paired FASTQ files.

ENA exposes paired compressed FASTQ files for eight runs, totaling approximately 66.65 GB. At the time of inspection, ENA did not provide usable read counts, base counts or FASTQ links for MS-1 (`SRR36739079`), although NCBI contains the corresponding run with 10.35 Gb of sequence. Extrapolation from the eight ENA FASTQ sets gives an estimated complete compressed FASTQ download of approximately **74 GB**. This is an estimate rather than a verified exact total. Storage required for decompressed FASTQ files, quality-control outputs and subsequent analysis products would be substantially greater, but no exact working-space requirement has yet been established.

The publication states that the primary data underlying its analyses are supplied as supplementary materials. Ten files are exposed through the article and PMC record: `DataSheet1.csv`, `DataSheet2.csv`, `Table1.xlsx` through `Table6.xlsx`, `Table7.xls` and `Table8.docx`. Their combined size is approximately 1.72 MB, with `Table4.xlsx` accounting for about 1.64 MB. The generic filenames do not by themselves establish the precise contents of every file, so their roles should remain provisional until the supplementary material is examined directly.

The BioProject record did not expose a separately deposited metagenome assembly, MAG collection, annotated gene catalogue or other dedicated processed-metagenomics repository. The identifiable processed or supporting data therefore consist of the publication supplements rather than a linked assembly or catalogue resource.

## 6. Metadata limitations and inconsistencies

The public metadata leave the following limitations unresolved.

- Treatment identity, salt composition and numerical salinity are not encoded as dedicated BioSample attributes. Treatment assignment is conveyed by the CK, LS and MS sample names and documented in the publication.
- BioSample does not contain per-sample soil chemistry or plant-trait measurements.
- Extraction batch, library batch, sequencing lane, experimental block, randomization procedure and pot identifiers are not reported in the repository metadata.
- The publication depicts eight pots per treatment but reports and deposits only three biological replicates per treatment. The selection of sequenced pots and the exact sample-to-pot relationship are unresolved.
- The material is operationally defined root-zone soil, not a strictly isolated rhizosphere fraction.
- Sampling depth and elevation were not collected.
- SRA RunInfo reports an insert size of zero even though the publication describes DNA fragmentation to approximately 450 bp.
- ENA currently lacks a usable FASTQ record for MS-1, whereas NCBI has the run and its sequence statistics.
- ENA’s `read_count` values for the available runs correspond to NCBI paired spots rather than to the doubled number of individual reads. Counts should therefore be labelled carefully.
- The BioProject-level archive size and the sum of run-level NCBI sizes differ. This likely reflects differing archive representations or record states, but the discrepancy has not been resolved.

## 7. Scope of the completed analyses

[Report 1](../report/alfalfa_salinity_metagenome_reanalysis.html) analyzes the nine source-derived sample profiles, their measured environmental context and 30 selected functional pathways. [Report 2](../report/alfalfa_ck1_raw_read_reconstruction.html) reconstructs CK-1 and MS-1 from NCBI SRA archives through QC, host depletion and Kraken2/Bracken taxonomic profiling. The two-sample comparison is descriptive; treatment-level inference remains limited to the supplementary-data analysis with three samples per treatment.

Neither report reconstructs the raw reads for all nine samples or a complete functional gene catalogue.

## 8. Source record

The metadata in this note were drawn from the following public records:

- [NCBI BioProject PRJNA1399466](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1399466)
- [NCBI BioSample records SAMN54495940–SAMN54495948](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=biosample&id=SAMN54495940,SAMN54495941,SAMN54495942,SAMN54495943,SAMN54495944,SAMN54495945,SAMN54495946,SAMN54495947,SAMN54495948&retmode=xml)
- [NCBI SRA RunInfo for PRJNA1399466](https://trace.ncbi.nlm.nih.gov/Traces/sra-db-be/runinfo?acc=PRJNA1399466)
- [ENA read-run file report for PRJNA1399466](https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA1399466&result=read_run&fields=run_accession,sample_accession,secondary_sample_accession,experiment_accession,study_accession,secondary_study_accession,sample_alias,library_name,library_strategy,library_source,library_selection,library_layout,instrument_platform,instrument_model,read_count,base_count,fastq_bytes,fastq_ftp&format=tsv)
- [Lv et al. (2026), *Frontiers in Plant Science*](https://www.frontiersin.org/journals/plant-science/articles/10.3389/fpls.2026.1753229/full), DOI [10.3389/fpls.2026.1753229](https://doi.org/10.3389/fpls.2026.1753229)
