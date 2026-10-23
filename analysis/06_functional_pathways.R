library(tidyverse)
library(here)
library(vegan)
source(here("analysis", "output_paths.R"), local = TRUE)

# ------------------------------------------------------------
# 1. Import validated pathway abundance matrix
# ------------------------------------------------------------

pathways <- read_tsv(
  here(
    "data",
    "processed",
    "pathway_abundance.tsv"
  ),
  show_col_types = FALSE
)

dim(pathways)
names(pathways)
pathways


# ------------------------------------------------------------
# 2. Basic inspection
# ------------------------------------------------------------

sample_cols <- names(pathways)[-1]

pathway_matrix <- as.matrix(
  pathways[sample_cols]
)

rownames(pathway_matrix) <- pathways$pathway

range(pathway_matrix)

sum(is.na(pathway_matrix))
sum(pathway_matrix < 0)
sum(pathway_matrix == 0)

sample_pathway_totals <- tibble(
  sample = sample_cols,
  total_reported_pathway_abundance = colSums(pathway_matrix)
)

sample_pathway_totals

# ------------------------------------------------------------
# 3. Pathway-level variability
# ------------------------------------------------------------

pathway_variability <- pathways |>
  rowwise() |>
  mutate(
    mean_abundance = mean(
      c_across(all_of(sample_cols))
    ),

    sd_abundance = sd(
      c_across(all_of(sample_cols))
    ),

    min_abundance = min(
      c_across(all_of(sample_cols))
    ),

    max_abundance = max(
      c_across(all_of(sample_cols))
    ),

    range_abundance = max_abundance - min_abundance,

    cv_percent = 100 * sd_abundance / mean_abundance
  ) |>
  ungroup() |>
  select(
    pathway,
    mean_abundance,
    sd_abundance,
    range_abundance,
    cv_percent
  ) |>
  arrange(desc(cv_percent))

pathway_variability

# ------------------------------------------------------------
# 4. Long-format pathway table
# ------------------------------------------------------------

pathway_long <- pathways |>
  pivot_longer(
    cols = all_of(sample_cols),
    names_to = "sample",
    values_to = "abundance"
  ) |>
  mutate(
    treatment = str_extract(
      sample,
      "^[A-Z]+"
    ),

    treatment = factor(
      treatment,
      levels = c(
        "CK",
        "LS",
        "MS"
      )
    )
  )

pathway_long

# ------------------------------------------------------------
# 5. Treatment-level pathway summaries
# ------------------------------------------------------------

pathway_summary <- pathway_long |>
  group_by(
    pathway,
    treatment
  ) |>
  summarise(
    mean_abundance = mean(abundance),
    sd_abundance = sd(abundance),
    .groups = "drop"
  )

pathway_summary


# ------------------------------------------------------------
# 6. Kruskal-Wallis tests for each pathway
# ------------------------------------------------------------

kw_results <- pathway_long |>
  group_by(pathway) |>
  group_modify(\(.x, .y) {
    test <- kruskal.test(
      abundance ~ treatment,
      data = .x
    )

    tibble(
      H = unname(test$statistic),
      df = unname(test$parameter),
      p_value = test$p.value
    )
  }) |>
  ungroup() |>
  mutate(
    p_BH = p.adjust(
      p_value,
      method = "BH"
    )
  ) |>
  arrange(p_value)

kw_results

pathway_treatment_means <- pathway_long |>
  group_by(pathway, treatment) |>
  summarise(
    mean_abundance = mean(abundance),
    median_abundance = median(abundance),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = treatment,
    values_from = c(
      mean_abundance,
      median_abundance
    )
  )

kw_results_full <- kw_results |>
  left_join(
    pathway_treatment_means,
    by = "pathway"
  )

kw_results_full


write_tsv(
  kw_results_full,
  result_path("tables", "pathway_kruskal_wallis_results.tsv")
)

kw_results_full <- kw_results_full |>
  mutate(
    epsilon_squared = pmax(
      0,
      (H - 3 + 1) / (9 - 3)
    )
  )

kw_results_full |>
  select(
    pathway,
    H,
    epsilon_squared,
    p_value,
    p_BH,
    mean_abundance_CK,
    mean_abundance_LS,
    mean_abundance_MS
  ) |>
  print(n = 30)


# ------------------------------------------------------------
# 7. Functional-profile community matrix
# ------------------------------------------------------------

functional_matrix <- pathways |>
  column_to_rownames("pathway") |>
  as.matrix() |>
  t()

dim(functional_matrix)
rownames(functional_matrix)
rowSums(functional_matrix)
functional_matrix[, 1:5]

# ------------------------------------------------------------
# 8. Bray-Curtis dissimilarity of functional profiles
# ------------------------------------------------------------

functional_bray <- vegdist(
  functional_matrix,
  method = "bray"
)

functional_bray_matrix <- as.matrix(functional_bray)

round(
  functional_bray_matrix,
  4
)

# ------------------------------------------------------------
# 9. PCoA of functional Bray-Curtis distances
# ------------------------------------------------------------

functional_pcoa <- wcmdscale(
  functional_bray,
  k = 2,
  eig = TRUE,
  add = FALSE
)

functional_eigenvalues <- functional_pcoa$eig

functional_eigenvalues

sum(functional_eigenvalues < 0)

functional_eigenvalues[
  functional_eigenvalues < 0
]

functional_positive_eig <-
  functional_eigenvalues[
    functional_eigenvalues > 0
  ]

functional_variance <- (functional_positive_eig /
  sum(functional_positive_eig)) *
  100

functional_variance[1:2]

functional_pcoa_scores <-
  as.data.frame(
    functional_pcoa$points
  )

colnames(
  functional_pcoa_scores
)[1:2] <- c(
  "PCoA1",
  "PCoA2"
)

functional_pcoa_scores <-
  functional_pcoa_scores |>
  rownames_to_column("sample") |>
  mutate(
    treatment = str_extract(
      sample,
      "^[A-Z]+"
    ),
    treatment = factor(
      treatment,
      levels = c(
        "CK",
        "LS",
        "MS"
      )
    )
  )

functional_pcoa_scores


# ------------------------------------------------------------
# 10. Metadata
# ------------------------------------------------------------

functional_metadata <- tibble(
  sample = rownames(
    functional_matrix
  ),
  treatment = str_extract(
    sample,
    "^[A-Z]+"
  )
) |>
  mutate(
    treatment = factor(
      treatment,
      levels = c(
        "CK",
        "LS",
        "MS"
      )
    )
  )

functional_metadata

identical(
  attr(
    functional_bray,
    "Labels"
  ),
  functional_metadata$sample
)


# ------------------------------------------------------------
# 11. PERMANOVA
# ------------------------------------------------------------

set.seed(1998)

functional_permanova <- adonis2(
  functional_bray ~ treatment,
  data = functional_metadata,
  permutations = 9999
)

functional_permanova

# ------------------------------------------------------------
# 12. PERMDISP
# ------------------------------------------------------------

functional_dispersion <- betadisper(
  functional_bray,
  group = functional_metadata$treatment,
  type = "median"
)

set.seed(1998)

functional_dispersion_test <- permutest(
  functional_dispersion,
  permutations = 9999
)

functional_dispersion_test

functional_dispersion_values <- tibble(
  sample = names(
    functional_dispersion$distances
  ),
  distance_to_median = functional_dispersion$distances
) |>
  mutate(
    treatment = str_extract(
      sample,
      "^[A-Z]+"
    )
  )

functional_dispersion_values

functional_dispersion_summary <-
  functional_dispersion_values |>
  group_by(treatment) |>
  summarise(
    mean_distance = mean(distance_to_median),

    sd_distance = sd(distance_to_median),

    .groups = "drop"
  )

functional_dispersion_summary

write_tsv(
  functional_pcoa_scores,
  result_path("tables", "functional_pcoa_scores.tsv")
)

write_tsv(
  tibble(
    axis = c(
      "PCoA1",
      "PCoA2"
    ),
    variance_explained = functional_variance[1:2]
  ),
  result_path("tables", "functional_pcoa_variance.tsv")
)

# ------------------------------------------------------------
# 13. Is PCoA1 associated with total selected-pathway abundance?
# ------------------------------------------------------------

functional_pcoa_check <- functional_pcoa_scores |>
  left_join(
    sample_pathway_totals,
    by = "sample"
  )

functional_pcoa_check |>
  select(
    sample,
    treatment,
    PCoA1,
    total_reported_pathway_abundance
  )

cor.test(
  functional_pcoa_check$PCoA1,
  functional_pcoa_check$total_reported_pathway_abundance,
  method = "pearson"
)


cor(
  functional_pcoa_check$PCoA1,
  functional_pcoa_check$total_reported_pathway_abundance
)

# ------------------------------------------------------------
# 14. Composition-only sensitivity analysis
# ------------------------------------------------------------

functional_composition <- sweep(
  functional_matrix,
  MARGIN = 1,
  STATS = rowSums(functional_matrix),
  FUN = "/"
)

rowSums(functional_composition)

functional_bray_composition <- vegdist(
  functional_composition,
  method = "bray"
)

functional_bray_composition_matrix <-
  as.matrix(functional_bray_composition)

round(
  functional_bray_composition_matrix,
  5
)

functional_pcoa_composition <- wcmdscale(
  functional_bray_composition,
  k = 2,
  eig = TRUE
)

eig_comp <- functional_pcoa_composition$eig

positive_eig_comp <- eig_comp[
  eig_comp > 0
]

variance_comp <- (positive_eig_comp /
  sum(positive_eig_comp)) *
  100

variance_comp[1:2]

set.seed(1998)

functional_permanova_composition <- adonis2(
  functional_bray_composition ~ treatment,
  data = functional_metadata,
  permutations = 9999
)

functional_permanova_composition

functional_dispersion_composition <- betadisper(
  functional_bray_composition,
  functional_metadata$treatment,
  type = "median"
)

set.seed(1998)

permutest(
  functional_dispersion_composition,
  permutations = 9999
)

# ------------------------------------------------------------
# 15. Publication-ready PCoA of normalized functional profiles
# ------------------------------------------------------------

# Extract sample coordinates
functional_comp_scores <- as.data.frame(
  functional_pcoa_composition$points
)

colnames(functional_comp_scores)[1:2] <- c(
  "PCoA1",
  "PCoA2"
)

functional_comp_scores <- functional_comp_scores |>
  rownames_to_column("sample") |>
  mutate(
    treatment = str_extract(
      sample,
      "^[A-Z]+"
    ),
    treatment = factor(
      treatment,
      levels = c("CK", "LS", "MS")
    )
  )

functional_comp_scores

# Treatment centroids
functional_comp_centroids <- functional_comp_scores |>
  group_by(treatment) |>
  summarise(
    PCoA1 = mean(PCoA1),
    PCoA2 = mean(PCoA2),
    .groups = "drop"
  )

functional_comp_centroids

treatment_cols <- c(
  "CK" = "#0072B2",
  "LS" = "#E69F00",
  "MS" = "#009E73"
)

p_functional_comp_pub <- ggplot(
  functional_comp_scores,
  aes(
    x = PCoA1,
    y = PCoA2,
    colour = treatment
  )
) +

  # Reference axes
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey72"
  ) +

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey72"
  ) +

  # Individual biological replicates
  geom_point(
    size = 4,
    alpha = 0.85
  ) +

  # Sample labels
  geom_text(
    aes(label = sample),
    nudge_y = 0.00035,
    size = 3.4,
    colour = "black",
    show.legend = FALSE
  ) +

  # Treatment centroids
  geom_point(
    data = functional_comp_centroids,
    aes(
      x = PCoA1,
      y = PCoA2,
      colour = treatment
    ),
    shape = 18,
    size = 5.2,
    show.legend = FALSE
  ) +

  scale_colour_manual(
    values = treatment_cols
  ) +

  labs(
    x = paste0(
      "PCoA1 (",
      round(variance_comp[1], 1),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(variance_comp[2], 1),
      "%)"
    ),
    colour = "Treatment"
  ) +

  # Preserve ordination geometry
  coord_equal() +

  theme_classic(base_size = 12) +

  theme(
    axis.text = element_text(
      colour = "black",
      size = 10
    ),

    axis.title = element_text(
      size = 12
    ),

    axis.title.x = element_text(
      margin = margin(t = 8)
    ),

    axis.title.y = element_text(
      margin = margin(r = 8)
    ),

    axis.line = element_line(
      linewidth = 0.5
    ),

    axis.ticks = element_line(
      linewidth = 0.4
    ),

    legend.position = "right",

    legend.title = element_text(
      size = 10
    ),

    legend.text = element_text(
      size = 10
    ),

    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )

p_functional_comp_pub


ggsave(
  result_path("figures", "normalized_functional_pathway_pcoa.png"),
  p_functional_comp_pub,
  width = 7.5,
  height = 6,
  units = "in",
  dpi = 600
)
