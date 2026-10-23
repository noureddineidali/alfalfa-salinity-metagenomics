library(tidyverse)
library(vegan)
library(here)
source(here("analysis", "output_paths.R"), local = TRUE)

# ------------------------------------------------------------
# 1. Import genus-level abundance matrix
# ------------------------------------------------------------

abund <- read_tsv(
  here("data", "processed", "genus_abundance.tsv"),
  show_col_types = FALSE
)

sample_cols <- names(abund)[-1]

# taxa x samples
abund_matrix <- as.matrix(abund[sample_cols])


# ------------------------------------------------------------
# 2. Remove features that are zero in every sample
# ------------------------------------------------------------

keep_taxa <- rowSums(abund_matrix) > 0

abund_matrix_filtered <- abund_matrix[keep_taxa, ]

dim(abund_matrix_filtered)


# ------------------------------------------------------------
# 3. Convert to within-sample relative abundance
# ------------------------------------------------------------

relative_matrix <- sweep(
  abund_matrix_filtered,
  MARGIN = 2,
  STATS = colSums(abund_matrix_filtered),
  FUN = "/"
)

colSums(relative_matrix)

# ------------------------------------------------------------
# 4. Samples as rows, taxa as columns
# ------------------------------------------------------------

community <- t(relative_matrix)

dim(community)
rownames(community)

# ------------------------------------------------------------
# 5. Bray-Curtis dissimilarity
# ------------------------------------------------------------

bray <- vegdist(
  community,
  method = "bray"
)

bray

# ------------------------------------------------------------
# 6. Inspect the full 9 x 9 Bray-Curtis matrix
# ------------------------------------------------------------

bray_matrix <- as.matrix(bray)

round(bray_matrix, 4)

# ------------------------------------------------------------
# 7. Bray-Curtis heatmap
# ------------------------------------------------------------

sample_order <- c(
  "CK-1",
  "CK-2",
  "CK-3",
  "LS-1",
  "LS-2",
  "LS-3",
  "MS-1",
  "MS-2",
  "MS-3"
)

# ------------------------------------------------------------
# 8. Unique pairwise comparisons
# ------------------------------------------------------------

distance_pairs <- as.data.frame(
  as.table(bray_matrix)
) |>
  rename(
    sample1 = Var1,
    sample2 = Var2,
    bray_curtis = Freq
  ) |>
  mutate(
    sample1 = as.character(sample1),
    sample2 = as.character(sample2)
  ) |>
  filter(
    match(sample1, sample_order) < match(sample2, sample_order)
  ) |>
  mutate(
    treatment1 = str_extract(sample1, "^[A-Z]+"),
    treatment2 = str_extract(sample2, "^[A-Z]+"),

    comparison = case_when(
      treatment1 == "CK" & treatment2 == "CK" ~ "Within CK",
      treatment1 == "LS" & treatment2 == "LS" ~ "Within LS",
      treatment1 == "MS" & treatment2 == "MS" ~ "Within MS",

      treatment1 %in% c("CK", "LS") & treatment2 %in% c("CK", "LS") ~ "CK-LS",

      treatment1 %in% c("CK", "MS") & treatment2 %in% c("CK", "MS") ~ "CK-MS",

      treatment1 %in% c("LS", "MS") & treatment2 %in% c("LS", "MS") ~ "LS-MS"
    ),

    comparison = factor(
      comparison,
      levels = c(
        "Within CK",
        "Within LS",
        "Within MS",
        "CK-LS",
        "CK-MS",
        "LS-MS"
      )
    )
  )


distance_pairs |>
  count(comparison)

distance_pairs |>
  group_by(comparison) |>
  summarise(
    mean = mean(bray_curtis),
    sd = sd(bray_curtis),
    .groups = "drop"
  )


# ------------------------------------------------------------
# Publication-ready pairwise Bray-Curtis comparison plot
# ------------------------------------------------------------

# ------------------------------------------------------------
# Mean ± SD for each comparison category
# ------------------------------------------------------------

distance_summary <- distance_pairs |>
  group_by(comparison) |>
  summarise(
    mean = mean(bray_curtis),
    sd = sd(bray_curtis),
    .groups = "drop"
  )


# ------------------------------------------------------------
# Colours
# ------------------------------------------------------------

comparison_cols <- c(
  "Within CK" = "#0072B2",
  "Within LS" = "#E69F00",
  "Within MS" = "#009E73",
  "CK-LS" = "#CC79A7",
  "CK-MS" = "#D55E00",
  "LS-MS" = "#56B4E9"
)


# ------------------------------------------------------------
# Plot
# ------------------------------------------------------------

p_bray_groups_pub <- ggplot(
  distance_pairs,
  aes(
    x = comparison,
    y = bray_curtis,
    colour = comparison
  )
) +

  # Individual pairwise Bray-Curtis distances
  geom_point(
    position = position_jitter(
      width = 0.08,
      height = 0,
      seed = 42
    ),
    size = 2.5,
    alpha = 0.75
  ) +

  # Mean ± SD
  geom_errorbar(
    data = distance_summary,
    aes(
      x = comparison,
      ymin = mean - sd,
      ymax = mean + sd,
      colour = comparison
    ),
    inherit.aes = FALSE,
    width = 0.14,
    linewidth = 0.65
  ) +

  # Mean
  geom_point(
    data = distance_summary,
    aes(
      x = comparison,
      y = mean,
      colour = comparison
    ),
    inherit.aes = FALSE,
    shape = 18,
    size = 4
  ) +

  # Separate within-treatment from between-treatment comparisons
  geom_vline(
    xintercept = 3.5,
    linetype = "dashed",
    colour = "grey60",
    linewidth = 0.5
  ) +

  scale_colour_manual(
    values = comparison_cols,
    guide = "none"
  ) +

  scale_y_continuous(
    breaks = seq(0.04, 0.18, by = 0.02),
    expand = expansion(mult = c(0.02, 0.04))
  ) +

  coord_cartesian(
    ylim = c(0.035, 0.195)
  ) +

  labs(
    x = NULL,
    y = "Bray-Curtis dissimilarity"
  ) +

  theme_classic(base_size = 12) +

  theme(
    axis.text.x = element_text(
      angle = 30,
      hjust = 1,
      colour = "black",
      size = 10
    ),

    axis.text.y = element_text(
      colour = "black",
      size = 10
    ),

    axis.title.y = element_text(
      size = 12,
      margin = margin(r = 8)
    ),

    axis.line = element_line(
      linewidth = 0.45
    ),

    axis.ticks = element_line(
      linewidth = 0.4
    ),

    plot.margin = margin(
      8,
      12,
      8,
      8
    )
  )

p_bray_groups_pub


# ------------------------------------------------------------
# Publication-ready triangular Bray-Curtis heatmap
# ------------------------------------------------------------

bray_heatmap <- as.data.frame(
  as.table(bray_matrix)
) |>
  rename(
    sample1 = Var1,
    sample2 = Var2,
    bray_curtis = Freq
  ) |>
  mutate(
    sample1 = as.character(sample1),
    sample2 = as.character(sample2),
    i = match(sample1, sample_order),
    j = match(sample2, sample_order)
  ) |>
  filter(i >= j) |>
  mutate(
    sample1 = factor(
      sample1,
      levels = rev(sample_order)
    ),
    sample2 = factor(
      sample2,
      levels = sample_order
    )
  )

p_bray_heatmap_pub <- ggplot(
  bray_heatmap,
  aes(
    x = sample2,
    y = sample1,
    fill = bray_curtis
  )
) +

  geom_tile(
    colour = "white",
    linewidth = 0.7
  ) +

  geom_text(
    aes(
      label = sprintf("%.3f", bray_curtis)
    ),
    size = 3.2,
    colour = "black"
  ) +

  scale_fill_gradient(
    low = "#F7FBFF",
    high = "#2166AC",
    limits = c(0, max(bray_matrix)),
    name = "Bray-Curtis\ndissimilarity"
  ) +

  # Treatment boundaries
  geom_vline(
    xintercept = c(3.5, 6.5),
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey45"
  ) +

  geom_hline(
    yintercept = c(3.5, 6.5),
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey45"
  ) +

  coord_fixed() +

  labs(
    x = NULL,
    y = NULL
  ) +

  theme_minimal(base_size = 11) +

  theme(
    panel.grid = element_blank(),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black",
      size = 10
    ),

    axis.text.y = element_text(
      colour = "black",
      size = 10
    ),

    axis.ticks = element_blank(),

    legend.title = element_text(
      size = 10
    ),

    legend.text = element_text(
      size = 9
    ),

    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )

p_bray_heatmap_pub


ggsave(
  result_path("figures", "bray_curtis_pairwise_groups.png"),
  p_bray_groups_pub,
  width = 8.5,
  height = 5.5,
  units = "in",
  dpi = 600
)


ggsave(
  result_path("figures", "bray_curtis_heatmap.png"),
  p_bray_heatmap_pub,
  width = 7.5,
  height = 6.5,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------
# 9. Principal Coordinates Analysis (PCoA)
# ------------------------------------------------------------

pcoa_fit <- wcmdscale(
  bray,
  k = 2,
  eig = TRUE,
  add = FALSE
)
pcoa_fit

eigenvalues <- pcoa_fit$eig
eigenvalues
sum(eigenvalues < 0)
eigenvalues[eigenvalues < 0]

variance_explained <- (eigenvalues /
  sum(eigenvalues)) *
  100

variance_explained[1:2]


# ------------------------------------------------------------
# 10. Extract PCoA coordinates
# ------------------------------------------------------------

pcoa_scores <- as.data.frame(pcoa_fit$points)

colnames(pcoa_scores)[1:2] <- c(
  "PCoA1",
  "PCoA2"
)

pcoa_scores <- pcoa_scores |>
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

pcoa_scores


# ------------------------------------------------------------
# 11. Publication-ready PCoA
# ------------------------------------------------------------

# Treatment centroids
pcoa_centroids <- pcoa_scores |>
  group_by(treatment) |>
  summarise(
    PCoA1 = mean(PCoA1),
    PCoA2 = mean(PCoA2),
    .groups = "drop"
  )

pcoa_centroids


treatment_cols <- c(
  "CK" = "#0072B2",
  "LS" = "#E69F00",
  "MS" = "#009E73"
)

p_pcoa_pub <- ggplot(
  pcoa_scores,
  aes(
    x = PCoA1,
    y = PCoA2,
    colour = treatment
  )
) +

  # Zero reference lines
  geom_hline(
    yintercept = 0,
    linewidth = 0.4,
    linetype = "dashed",
    colour = "grey70"
  ) +

  geom_vline(
    xintercept = 0,
    linewidth = 0.4,
    linetype = "dashed",
    colour = "grey70"
  ) +

  # Individual biological samples
  geom_point(
    size = 4,
    alpha = 0.85
  ) +

  # Sample labels
  geom_text(
    aes(label = sample),
    nudge_y = 0.006,
    size = 3.5,
    show.legend = FALSE
  ) +

  # Treatment centroids
  geom_point(
    data = pcoa_centroids,
    aes(
      x = PCoA1,
      y = PCoA2,
      colour = treatment
    ),
    shape = 18,
    size = 5,
    show.legend = FALSE
  ) +

  scale_colour_manual(
    values = treatment_cols
  ) +

  labs(
    x = paste0(
      "PCoA1 (",
      round(variance_explained[1], 1),
      "%)"
    ),
    y = paste0(
      "PCoA2 (",
      round(variance_explained[2], 1),
      "%)"
    ),
    colour = "Treatment"
  ) +

  # Important for ordination geometry
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

p_pcoa_pub


ggsave(
  result_path("figures", "genus_bray_curtis_pcoa.png"),
  p_pcoa_pub,
  width = 7.5,
  height = 6,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------
# 12. Sample metadata for multivariate tests
# ------------------------------------------------------------

metadata <- tibble(
  sample = rownames(community),
  treatment = str_extract(sample, "^[A-Z]+")
) |>
  mutate(
    treatment = factor(
      treatment,
      levels = c("CK", "LS", "MS")
    )
  )

metadata

# ------------------------------------------------------------
# 13. PERMANOVA
# ------------------------------------------------------------

set.seed(1998)

permanova <- adonis2(
  bray ~ treatment,
  data = metadata,
  permutations = 9999,
  by = "margin"
)

permanova


# ------------------------------------------------------------
# 14. PERMDISP
# ------------------------------------------------------------

dispersion <- betadisper(
  bray,
  group = metadata$treatment,
  type = "median"
)

dispersion

set.seed(1998)

dispersion_test <- permutest(
  dispersion,
  permutations = 9999
)

dispersion_test

# ------------------------------------------------------------
# 15. Mean distance to treatment median
# ------------------------------------------------------------

dispersion_values <- tibble(
  sample = names(dispersion$distances),
  distance_to_median = dispersion$distances
) |>
  mutate(
    treatment = str_extract(sample, "^[A-Z]+")
  )

dispersion_values


dispersion_summary <- dispersion_values |>
  group_by(treatment) |>
  summarise(
    mean_distance = mean(distance_to_median),
    sd_distance = sd(distance_to_median),
    .groups = "drop"
  )

dispersion_summary


# ------------------------------------------------------------
# Pairwise PERMANOVA: CK vs LS
# ------------------------------------------------------------

meta_ck_ls <- metadata |>
  filter(treatment %in% c("CK", "LS")) |>
  droplevels()

bray_ck_ls <- as.dist(
  bray_matrix[
    meta_ck_ls$sample,
    meta_ck_ls$sample
  ]
)

set.seed(1998)

perm_ck_ls <- adonis2(
  bray_ck_ls ~ treatment,
  data = meta_ck_ls,
  permutations = 9999
)

perm_ck_ls


# ------------------------------------------------------------
# Pairwise PERMANOVA: CK vs MS
# ------------------------------------------------------------

meta_ck_ms <- metadata |>
  filter(treatment %in% c("CK", "MS")) |>
  droplevels()

bray_ck_ms <- as.dist(
  bray_matrix[
    meta_ck_ms$sample,
    meta_ck_ms$sample
  ]
)

set.seed(1998)

perm_ck_ms <- adonis2(
  bray_ck_ms ~ treatment,
  data = meta_ck_ms,
  permutations = 9999
)

perm_ck_ms

# ------------------------------------------------------------
# Pairwise PERMANOVA: LS vs MS
# ------------------------------------------------------------

meta_ls_ms <- metadata |>
  filter(treatment %in% c("LS", "MS")) |>
  droplevels()

bray_ls_ms <- as.dist(
  bray_matrix[
    meta_ls_ms$sample,
    meta_ls_ms$sample
  ]
)

set.seed(1998)

perm_ls_ms <- adonis2(
  bray_ls_ms ~ treatment,
  data = meta_ls_ms,
  permutations = 9999
)

perm_ls_ms

pairwise_permanova <- tibble(
  comparison = c("CK-LS", "CK-MS", "LS-MS"),

  R2 = c(
    perm_ck_ls$R2[1],
    perm_ck_ms$R2[1],
    perm_ls_ms$R2[1]
  ),

  pseudo_F = c(
    perm_ck_ls$F[1],
    perm_ck_ms$F[1],
    perm_ls_ms$F[1]
  ),

  p_value = c(
    perm_ck_ls$`Pr(>F)`[1],
    perm_ck_ms$`Pr(>F)`[1],
    perm_ls_ms$`Pr(>F)`[1]
  )
)

pairwise_permanova <- pairwise_permanova |>
  mutate(
    p_BH = p.adjust(p_value, method = "BH")
  )

pairwise_permanova

# ------------------------------------------------------------
# Export PCoA results for downstream analyses
# ------------------------------------------------------------

write_tsv(
  pcoa_scores,
  result_path("tables", "pcoa_scores.tsv")
)

write_tsv(
  tibble(
    axis = c("PCoA1", "PCoA2"),
    variance_explained = variance_explained[1:2]
  ),
  result_path("tables", "pcoa_variance_explained.tsv")
)
