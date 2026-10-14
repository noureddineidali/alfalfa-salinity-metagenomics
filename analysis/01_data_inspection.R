library(tidyverse)
library(here)

# ------------------------------------------------------------
# 1. Import the published genus-level abundance matrix
# ------------------------------------------------------------

abund <- read_tsv(
  here("data", "processed", "genus_abundance.tsv"),
  show_col_types = FALSE
)

# Basic structure
dim(abund)
names(abund)
glimpse(abund)


# ------------------------------------------------------------
# 2. Identify taxon and sample columns
# ------------------------------------------------------------

taxon_col <- names(abund)[1]
sample_cols <- names(abund)[-1]

taxon_col
sample_cols

# We expect exactly 9 samples
stopifnot(length(sample_cols) == 9)


# ------------------------------------------------------------
# 3. Basic data-quality checks
# ------------------------------------------------------------

# Duplicated reported taxon labels
sum(duplicated(abund[[taxon_col]]))

abund |>
  distinct(genus) |>
  count()

# Missing values
sum(is.na(abund))

# Negative abundance values
sum(as.matrix(abund[sample_cols]) < 0, na.rm = TRUE)


# ------------------------------------------------------------
# 4. Sample-level summary
# ------------------------------------------------------------

abund_matrix <- as.matrix(abund[sample_cols])

sample_summary <- tibble(
  sample = sample_cols,

  # Sum of the undocumented supplementary abundance values
  total_abundance = colSums(abund_matrix, na.rm = TRUE),

  # Number of reported genus-level labels with abundance > 0
  observed_features = colSums(abund_matrix > 0, na.rm = TRUE),

  # Number of zero-valued features
  zero_features = colSums(abund_matrix == 0, na.rm = TRUE),

  # Missing values
  missing_values = colSums(is.na(abund_matrix)),

  # Negative values
  negative_values = colSums(abund_matrix < 0, na.rm = TRUE)
) |>
  mutate(
    zero_fraction = zero_features / nrow(abund)
  )

sample_summary


# ------------------------------------------------------------
# 5. How sparse is the matrix?
# ------------------------------------------------------------

overall_zero_fraction <-
  sum(abund_matrix == 0, na.rm = TRUE) /
  length(abund_matrix)

overall_zero_fraction


# ------------------------------------------------------------
# 6. Prevalence of each reported taxon label
# ------------------------------------------------------------

# In how many of the 9 samples is each feature detected (>0)?
feature_prevalence <- rowSums(abund_matrix > 0, na.rm = TRUE)

table(feature_prevalence)


# ------------------------------------------------------------
# 7. Convert abundance values to within-sample proportions
# ------------------------------------------------------------

relative_matrix <- sweep(
  abund_matrix,
  MARGIN = 2,
  STATS = colSums(abund_matrix, na.rm = TRUE),
  FUN = "/"
)

# Each sample should now sum to approximately 1
colSums(relative_matrix)


# Put the taxon labels back
relative_abund <- bind_cols(
  abund[1],
  as_tibble(relative_matrix)
)

names(relative_abund) <- names(abund)

# Inspect
relative_abund


# ------------------------------------------------------------
# 8. Summarize matrix structure by treatment
# ------------------------------------------------------------

sample_summary_by_treatment <- sample_summary |>
  mutate(
    treatment = str_extract(sample, "^[A-Z]+")
  ) |>
  group_by(treatment) |>
  summarise(
    mean_observed_features = mean(observed_features),
    sd_observed_features = sd(observed_features),
    mean_zero_fraction = mean(zero_fraction),
    sd_zero_fraction = sd(zero_fraction),
    .groups = "drop"
  )

sample_summary_by_treatment


# ------------------------------------------------------------
# 9. Convert relative abundance table to long format
# ------------------------------------------------------------

relative_long <- relative_abund |>
  rename(taxon = all_of(taxon_col)) |>
  pivot_longer(
    cols = -taxon,
    names_to = "sample",
    values_to = "relative_abundance"
  ) |>
  mutate(
    treatment = str_extract(sample, "^[A-Z]+"),
    replicate = as.integer(str_extract(sample, "\\d+$"))
  )

relative_long

# ------------------------------------------------------------
# 10. Taxon prevalence and mean abundance within treatments
# ------------------------------------------------------------

treatment_taxon_summary <- relative_long |>
  group_by(taxon, treatment) |>
  summarise(
    prevalence_n = sum(relative_abundance > 0),
    prevalence_prop = mean(relative_abundance > 0),
    mean_relative_abundance = mean(relative_abundance),
    median_relative_abundance = median(relative_abundance),
    .groups = "drop"
  )

treatment_taxon_summary

# ------------------------------------------------------------
# 11. Compare prevalence patterns across CK, LS and MS
# ------------------------------------------------------------

prevalence_wide <- treatment_taxon_summary |>
  select(taxon, treatment, prevalence_n) |>
  pivot_wider(
    names_from = treatment,
    values_from = prevalence_n
  )

prevalence_wide

prevalence_patterns <- prevalence_wide |>
  count(CK, LS, MS, sort = TRUE)

prevalence_patterns
