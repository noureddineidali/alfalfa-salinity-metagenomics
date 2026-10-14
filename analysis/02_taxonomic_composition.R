library(tidyverse)
library(here)
library(scales)

# ------------------------------------------------------------
# 1. Import abundance matrix
# ------------------------------------------------------------

abund <- read_tsv(
  here("data", "processed", "genus_abundance.tsv"),
  show_col_types = FALSE
)

taxon_col <- names(abund)[1]
sample_cols <- names(abund)[-1]

abund_matrix <- as.matrix(abund[sample_cols])

# ------------------------------------------------------------
# 2. Convert to relative abundance
# ------------------------------------------------------------

relative_matrix <- sweep(
  abund_matrix,
  MARGIN = 2,
  STATS = colSums(abund_matrix),
  FUN = "/"
)

relative_abund <- bind_cols(
  taxon = abund[[taxon_col]],
  as_tibble(relative_matrix)
)

# ------------------------------------------------------------
# 3. Long format
# ------------------------------------------------------------

relative_long <- relative_abund |>
  pivot_longer(
    cols = -taxon,
    names_to = "sample",
    values_to = "relative_abundance"
  ) |>
  mutate(
    treatment = str_extract(sample, "^[A-Z]+"),
    replicate = as.integer(str_extract(sample, "\\d+$"))
  )

# ------------------------------------------------------------
# 4. Rank taxa by mean relative abundance
# ------------------------------------------------------------

taxon_ranking <- relative_long |>
  group_by(taxon) |>
  summarise(
    mean_relative_abundance = mean(relative_abundance),
    .groups = "drop"
  ) |>
  arrange(desc(mean_relative_abundance))

top10_taxa <- taxon_ranking |>
  slice_head(n = 10) |>
  pull(taxon)

# ------------------------------------------------------------
# 5. Collapse everything else into Others
# ------------------------------------------------------------

composition_plot_data <- relative_long |>
  mutate(
    display_taxon = if_else(
      taxon %in% top10_taxa,
      taxon,
      "Others"
    )
  ) |>
  group_by(sample, treatment, display_taxon) |>
  summarise(
    relative_abundance = sum(relative_abundance),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 6. Sample order
# ------------------------------------------------------------

sample_levels <- c(
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

composition_plot_data <- composition_plot_data |>
  mutate(
    sample = factor(sample, levels = sample_levels)
  )

# ------------------------------------------------------------
# 7. Taxon order for plotting
# ------------------------------------------------------------

taxon_levels <- c(rev(top10_taxa), "Others")

composition_plot_data <- composition_plot_data |>
  mutate(
    display_taxon = factor(display_taxon, levels = taxon_levels)
  )

# ------------------------------------------------------------
# 8. Clean labels if desired
# ------------------------------------------------------------

label_map <- c(
  "SCGC-AG-212-J23" = "SCGC-AG-212-J23",
  "Gp6-AA40" = "Gp6-AA40",
  "Sediminibacterium" = "Sediminibacterium",
  "PSRF01" = "PSRF01",
  "12-FULL-67-14b" = "12-FULL-67-14b",
  "Bradyrhizobium" = "Bradyrhizobium",
  "UBA12294" = "UBA12294",
  "JACDDX01" = "JACDDX01",
  "WHTF01" = "WHTF01",
  "Nocardioides" = "Nocardioides",
  "Others" = "Others"
)

# ------------------------------------------------------------
# 9. Publication-style composition plot
# ------------------------------------------------------------

p1_pub <- ggplot(
  composition_plot_data,
  aes(
    x = sample,
    y = relative_abundance,
    fill = display_taxon
  )
) +
  geom_col(width = 0.82, colour = "white", linewidth = 0.15) +
  geom_vline(
    xintercept = c(3.5, 6.5),
    linewidth = 0.4,
    linetype = "dashed",
    colour = "grey40"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_fill_brewer(
    palette = "Paired",
    breaks = taxon_levels,
    labels = label_map
  ) +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Reported genus-level taxon"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(colour = "black"),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.position = "right",
    panel.border = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

p1_pub


# ggsave(
#   filename = here(
#     "results",
#     "figures",
#     "top10_plus_others_by_sample.pdf"
#   ),
#   plot = p1_pub,
#   width = 9,
#   height = 6,
#   units = "in"
# )

# ggsave(
#   filename = here(
#     "results",
#     "figures",
#     "top10_plus_others_by_sample.png"
#   ),
#   plot = p1_pub,
#   width = 9,
#   height = 6,
#   units = "in",
#   dpi = 600
# )

# ------------------------------------------------------------
# 10. Top-10 taxa only — true whole-community abundance
# ------------------------------------------------------------

top10_plot_data <- relative_long |>
  filter(taxon %in% top10_taxa) |>
  mutate(
    sample = factor(
      sample,
      levels = c(
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
    ),
    taxon = factor(
      taxon,
      levels = rev(top10_taxa)
    )
  )

p2_top10 <- ggplot(
  top10_plot_data,
  aes(
    x = sample,
    y = relative_abundance,
    fill = taxon
  )
) +
  geom_col(
    width = 0.82,
    colour = "white",
    linewidth = 0.15
  ) +
  geom_vline(
    xintercept = c(3.5, 6.5),
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_fill_brewer(
    palette = "Paired"
  ) +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Reported genus-level taxon"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(colour = "black"),
    legend.position = "right"
  )

p2_top10


# ggsave(
#   here(
#     "results",
#     "figures",
#     "top10_taxa_by_sample.pdf"
#   ),
#   p2_top10,
#   width = 9,
#   height = 6,
#   units = "in"
# )

# ggsave(
#   here(
#     "results",
#     "figures",
#     "top10_taxa_by_sample.png"
#   ),
#   p2_top10,
#   width = 9,
#   height = 6,
#   units = "in",
#   dpi = 600
# )

# ------------------------------------------------------------
# 11. Treatment-level summary of top-10 taxa
# ------------------------------------------------------------

top10_treatment_summary <- relative_long |>
  filter(taxon %in% top10_taxa) |>
  group_by(taxon, treatment) |>
  summarise(
    mean_relative_abundance = mean(relative_abundance),
    sd_relative_abundance = sd(relative_abundance),
    median_relative_abundance = median(relative_abundance),
    min_relative_abundance = min(relative_abundance),
    max_relative_abundance = max(relative_abundance),
    .groups = "drop"
  )

top10_treatment_summary

top10_treatment_summary_percent <- top10_treatment_summary |>
  mutate(
    across(
      c(
        mean_relative_abundance,
        sd_relative_abundance,
        median_relative_abundance,
        min_relative_abundance,
        max_relative_abundance
      ),
      ~ .x * 100
    )
  )

top10_treatment_summary_percent

# ------------------------------------------------------------
# 12. Individual replicates + treatment mean ± SD
# ------------------------------------------------------------

library(tidyverse)
library(here)
library(scales)

# ------------------------------------------------------------
# 1. Import abundance matrix
# ------------------------------------------------------------

abund <- read_tsv(
  here("data", "processed", "genus_abundance.tsv"),
  show_col_types = FALSE
)

taxon_col <- names(abund)[1]
sample_cols <- names(abund)[-1]

abund_matrix <- as.matrix(abund[sample_cols])

# ------------------------------------------------------------
# 2. Convert to relative abundance
# ------------------------------------------------------------

relative_matrix <- sweep(
  abund_matrix,
  MARGIN = 2,
  STATS = colSums(abund_matrix),
  FUN = "/"
)

relative_abund <- bind_cols(
  taxon = abund[[taxon_col]],
  as_tibble(relative_matrix)
)

# ------------------------------------------------------------
# 3. Long format
# ------------------------------------------------------------

relative_long <- relative_abund |>
  pivot_longer(
    cols = -taxon,
    names_to = "sample",
    values_to = "relative_abundance"
  ) |>
  mutate(
    treatment = str_extract(sample, "^[A-Z]+"),
    replicate = as.integer(str_extract(sample, "\\d+$"))
  )

# ------------------------------------------------------------
# 4. Rank taxa by mean relative abundance
# ------------------------------------------------------------

taxon_ranking <- relative_long |>
  group_by(taxon) |>
  summarise(
    mean_relative_abundance = mean(relative_abundance),
    .groups = "drop"
  ) |>
  arrange(desc(mean_relative_abundance))

top10_taxa <- taxon_ranking |>
  slice_head(n = 10) |>
  pull(taxon)

# ------------------------------------------------------------
# 5. Collapse everything else into Others
# ------------------------------------------------------------

composition_plot_data <- relative_long |>
  mutate(
    display_taxon = if_else(
      taxon %in% top10_taxa,
      taxon,
      "Others"
    )
  ) |>
  group_by(sample, treatment, display_taxon) |>
  summarise(
    relative_abundance = sum(relative_abundance),
    .groups = "drop"
  )

# ------------------------------------------------------------
# 6. Sample order
# ------------------------------------------------------------

sample_levels <- c(
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

composition_plot_data <- composition_plot_data |>
  mutate(
    sample = factor(sample, levels = sample_levels)
  )

# ------------------------------------------------------------
# 7. Taxon order for plotting
# ------------------------------------------------------------

taxon_levels <- c(rev(top10_taxa), "Others")

composition_plot_data <- composition_plot_data |>
  mutate(
    display_taxon = factor(display_taxon, levels = taxon_levels)
  )

# ------------------------------------------------------------
# 8. Clean labels if desired
# ------------------------------------------------------------

label_map <- c(
  "SCGC-AG-212-J23" = "SCGC-AG-212-J23",
  "Gp6-AA40" = "Gp6-AA40",
  "Sediminibacterium" = "Sediminibacterium",
  "PSRF01" = "PSRF01",
  "12-FULL-67-14b" = "12-FULL-67-14b",
  "Bradyrhizobium" = "Bradyrhizobium",
  "UBA12294" = "UBA12294",
  "JACDDX01" = "JACDDX01",
  "WHTF01" = "WHTF01",
  "Nocardioides" = "Nocardioides",
  "Others" = "Others"
)

# ------------------------------------------------------------
# 9. Publication-style composition plot
# ------------------------------------------------------------

p1_pub <- ggplot(
  composition_plot_data,
  aes(
    x = sample,
    y = relative_abundance,
    fill = display_taxon
  )
) +
  geom_col(width = 0.82, colour = "white", linewidth = 0.15) +
  geom_vline(
    xintercept = c(3.5, 6.5),
    linewidth = 0.4,
    linetype = "dashed",
    colour = "grey40"
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.01))
  ) +
  scale_fill_brewer(
    palette = "Paired",
    breaks = taxon_levels,
    labels = label_map
  ) +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Reported genus-level taxon"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(colour = "black"),
    axis.title.y = element_text(size = 12),
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 9),
    legend.position = "right",
    panel.border = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

p1_pub


# ggsave(
#   filename = here(
#     "results",
#     "figures",
#     "top10_plus_others_by_sample.pdf"
#   ),
#   plot = p1_pub,
#   width = 9,
#   height = 6,
#   units = "in"
# )

# ggsave(
#   filename = here(
#     "results",
#     "figures",
#     "top10_plus_others_by_sample.png"
#   ),
#   plot = p1_pub,
#   width = 9,
#   height = 6,
#   units = "in",
#   dpi = 600
# )

# ------------------------------------------------------------
# 10. Top-10 taxa only — true whole-community abundance
# ------------------------------------------------------------

top10_plot_data <- relative_long |>
  filter(taxon %in% top10_taxa) |>
  mutate(
    sample = factor(
      sample,
      levels = c(
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
    ),
    taxon = factor(
      taxon,
      levels = rev(top10_taxa)
    )
  )

p2_top10 <- ggplot(
  top10_plot_data,
  aes(
    x = sample,
    y = relative_abundance,
    fill = taxon
  )
) +
  geom_col(
    width = 0.82,
    colour = "white",
    linewidth = 0.15
  ) +
  geom_vline(
    xintercept = c(3.5, 6.5),
    linetype = "dashed",
    linewidth = 0.4,
    colour = "grey40"
  ) +
  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    expand = expansion(mult = c(0, 0.03))
  ) +
  scale_fill_brewer(
    palette = "Paired"
  ) +
  labs(
    x = NULL,
    y = "Relative abundance",
    fill = "Reported genus-level taxon"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      colour = "black"
    ),
    axis.text.y = element_text(colour = "black"),
    legend.position = "right"
  )

p2_top10


# ggsave(
#   here(
#     "results",
#     "figures",
#     "top10_taxa_by_sample.pdf"
#   ),
#   p2_top10,
#   width = 9,
#   height = 6,
#   units = "in"
# )

# ggsave(
#   here(
#     "results",
#     "figures",
#     "top10_taxa_by_sample.png"
#   ),
#   p2_top10,
#   width = 9,
#   height = 6,
#   units = "in",
#   dpi = 600
# )

# ------------------------------------------------------------
# 11. Treatment-level summary of top-10 taxa
# ------------------------------------------------------------

top10_treatment_summary <- relative_long |>
  filter(taxon %in% top10_taxa) |>
  group_by(taxon, treatment) |>
  summarise(
    mean_relative_abundance = mean(relative_abundance),
    sd_relative_abundance = sd(relative_abundance),
    median_relative_abundance = median(relative_abundance),
    min_relative_abundance = min(relative_abundance),
    max_relative_abundance = max(relative_abundance),
    .groups = "drop"
  )

top10_treatment_summary

top10_treatment_summary_percent <- top10_treatment_summary |>
  mutate(
    across(
      c(
        mean_relative_abundance,
        sd_relative_abundance,
        median_relative_abundance,
        min_relative_abundance,
        max_relative_abundance
      ),
      ~ .x * 100
    )
  )

top10_treatment_summary_percent

# ------------------------------------------------------------
# 12. Individual replicates + treatment mean ± SD
# ------------------------------------------------------------

top10_points <- relative_long |>
  filter(taxon %in% top10_taxa) |>
  mutate(
    treatment = factor(
      treatment,
      levels = c("CK", "LS", "MS")
    ),
    taxon = factor(
      taxon,
      levels = top10_taxa
    )
  )


# Explicit summary statistics
top10_summary <- top10_points |>
  group_by(taxon, treatment) |>
  summarise(
    mean_abundance = mean(relative_abundance),
    sd_abundance = sd(relative_abundance),
    .groups = "drop"
  ) |>
  mutate(
    ymin = pmax(mean_abundance - sd_abundance, 0),
    ymax = mean_abundance + sd_abundance
  )


# Colourblind-safe treatment palette
treatment_cols <- c(
  "CK" = "#0072B2",
  "LS" = "#E69F00",
  "MS" = "#009E73"
)

p3_treatment <- ggplot(
  top10_points,
  aes(
    x = treatment,
    y = relative_abundance
  )
) +

  # Individual biological replicates
  geom_point(
    aes(colour = treatment),
    position = position_jitter(
      width = 0.06,
      height = 0,
      seed = 42
    ),
    size = 2.2,
    alpha = 0.70
  ) +

  # Mean ± SD
  geom_errorbar(
    data = top10_summary,
    aes(
      x = treatment,
      ymin = ymin,
      ymax = ymax,
      colour = treatment
    ),
    inherit.aes = FALSE,
    width = 0.13,
    linewidth = 0.65
  ) +

  geom_point(
    data = top10_summary,
    aes(
      x = treatment,
      y = mean_abundance,
      colour = treatment
    ),
    inherit.aes = FALSE,
    shape = 18,
    size = 3.5
  ) +

  facet_wrap(
    ~taxon,
    ncol = 5
  ) +

  scale_colour_manual(
    values = treatment_cols,
    guide = "none"
  ) +

  scale_y_continuous(
    labels = scales::percent_format(accuracy = 1),
    limits = c(0, NA),
    expand = expansion(mult = c(0, 0.06))
  ) +

  labs(
    x = NULL,
    y = "Relative abundance"
  ) +
  theme_classic(base_size = 11)

p3_treatment

# ggsave(
#   here(
#     "results",
#     "figures",
#     "top10_taxa_treatment_patterns.pdf"
#   ),
#   p3_treatment,
#   width = 11,
#   height = 6.5,
#   units = "in"
# )

# ggsave(
#   here(
#     "results",
#     "figures",
#     "top10_taxa_treatment_patterns.png"
#   ),
#   p3_treatment,
#   width = 11,
#   height = 6.5,
#   units = "in",
#   dpi = 600
# )
