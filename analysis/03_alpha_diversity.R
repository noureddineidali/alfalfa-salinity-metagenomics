library(tidyverse)
library(vegan)
library(here)

# ------------------------------------------------------------
# 1. Import abundance matrix
# ------------------------------------------------------------

abund <- read_tsv(
  here("data", "processed", "genus_abundance.tsv"),
  show_col_types = FALSE
)

sample_cols <- names(abund)[-1]

# vegan expects samples as rows and taxa as columns
community <- t(
  as.matrix(abund[sample_cols])
)

# ------------------------------------------------------------
# 2. Alpha-diversity metrics
# ------------------------------------------------------------

alpha_diversity <- tibble(
  sample = rownames(community),

  richness = specnumber(community),

  shannon = diversity(
    community,
    index = "shannon"
  ),

  simpson = diversity(
    community,
    index = "simpson"
  ),

  inverse_simpson = diversity(
    community,
    index = "invsimpson"
  )
) |>
  mutate(
    treatment = str_extract(sample, "^[A-Z]+"),
    treatment = factor(
      treatment,
      levels = c("CK", "LS", "MS")
    ),

    # Convert Shannon from natural-log base to log2
    shannon_log2 = shannon / log(2)
  )

alpha_diversity


# ------------------------------------------------------------
# 3. Prepare publication plot data
# ------------------------------------------------------------

alpha_plot_data <- alpha_diversity |>
  select(
    sample,
    treatment,
    richness,
    shannon_log2,
    simpson
  ) |>
  pivot_longer(
    cols = c(
      richness,
      shannon_log2,
      simpson
    ),
    names_to = "metric",
    values_to = "value"
  ) |>
  mutate(
    metric = factor(
      metric,
      levels = c(
        "richness",
        "shannon_log2",
        "simpson"
      ),
      labels = c(
        "Observed richness",
        "Shannon diversity",
        "Simpson diversity"
      )
    )
  )


# ------------------------------------------------------------
# 4. Explicit means and SD
# ------------------------------------------------------------

alpha_plot_summary <- alpha_plot_data |>
  group_by(metric, treatment) |>
  summarise(
    mean = mean(value),
    sd = sd(value),
    .groups = "drop"
  )


# ------------------------------------------------------------
# 5. Treatment colours
# ------------------------------------------------------------

treatment_cols <- c(
  "CK" = "#0072B2",
  "LS" = "#E69F00",
  "MS" = "#009E73"
)


# ------------------------------------------------------------
# 6. Publication-ready alpha-diversity figure
# ------------------------------------------------------------

p_alpha_pub <- ggplot(
  alpha_plot_data,
  aes(
    x = treatment,
    y = value,
    colour = treatment
  )
) +

  # Individual biological replicates
  geom_point(
    position = position_jitter(
      width = 0.06,
      height = 0,
      seed = 42
    ),
    size = 2.5,
    alpha = 0.75
  ) +

  # Mean ± SD
  geom_errorbar(
    data = alpha_plot_summary,
    aes(
      x = treatment,
      ymin = mean - sd,
      ymax = mean + sd,
      colour = treatment
    ),
    inherit.aes = FALSE,
    width = 0.14,
    linewidth = 0.65
  ) +

  # Mean
  geom_point(
    data = alpha_plot_summary,
    aes(
      x = treatment,
      y = mean,
      colour = treatment
    ),
    inherit.aes = FALSE,
    shape = 18,
    size = 3.7
  ) +

  facet_wrap(
    ~metric,
    nrow = 1,
    scales = "free_y"
  ) +

  scale_colour_manual(
    values = treatment_cols,
    guide = "none"
  ) +

  labs(
    x = NULL,
    y = NULL
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    # Clean facet titles
    strip.background = element_blank(),

    strip.text = element_text(
      face = "bold",
      size = 11,
      colour = "black",
      margin = margin(
        b = 7
      )
    ),

    axis.text.x = element_text(
      colour = "black",
      size = 10
    ),

    axis.text.y = element_text(
      colour = "black",
      size = 9
    ),

    axis.line = element_line(
      linewidth = 0.45
    ),

    axis.ticks = element_line(
      linewidth = 0.4
    ),

    panel.spacing = unit(
      1.4,
      "lines"
    ),

    plot.margin = margin(
      8,
      8,
      8,
      8
    )
  )

p_alpha_pub

# ------------------------------------------------------------
# 7. Save figure
# ------------------------------------------------------------

# ggsave(
#   filename = here(
#     "results",
#     "figures",
#     "alpha_diversity_genus_level.pdf"
#   ),
#   plot = p_alpha_pub,
#   width = 10,
#   height = 4.5,
#   units = "in"
# )

# ggsave(
#   filename = here(
#     "results",
#     "figures",
#     "alpha_diversity_genus_level.png"
#   ),
#   plot = p_alpha_pub,
#   width = 10,
#   height = 4.5,
#   units = "in",
#   dpi = 600
# )

# ------------------------------------------------------------
# Save genus-level alpha-diversity results
# ------------------------------------------------------------

alpha_results <- alpha_diversity |>
  mutate(
    shannon_log2 = shannon / log(2)
  ) |>
  select(
    sample,
    treatment,
    richness,
    shannon_ln = shannon,
    shannon_log2,
    simpson,
    inverse_simpson
  )

# write_tsv(
#   alpha_results,
#   here(
#     "results",
#     "tables",
#     "genus_alpha_diversity.tsv"
#   )
# )

alpha_results
