library(tidyverse)
library(here)
library(vegan)
library(janitor)
source(here("analysis", "output_paths.R"), local = TRUE)

# ------------------------------------------------------------
# 1. Import environmental + plant metadata
# ------------------------------------------------------------

env <- read_tsv(
  here(
    "data",
    "processed",
    "environmental_plant_metadata.tsv"
  ),
  show_col_types = FALSE
) |>
  clean_names()

glimpse(env)
env

env <- env |>
  mutate(
    treatment_code = case_when(
      treatment == "Control" ~ "CK",
      treatment == "Low salinity" ~ "LS",
      treatment == "Moderate salinity" ~ "MS"
    ),
    treatment_code = factor(
      treatment_code,
      levels = c("CK", "LS", "MS")
    )
  )

# ------------------------------------------------------------
# 2. Prepare soil variables
# ------------------------------------------------------------

soil_vars <- c(
  "ec_us_cm",
  "ph",
  "total_n_g_kg",
  "nitrate_n_mg_kg",
  "ammonium_n_mg_kg"
)

# Standardize each soil variable
soil_z <- env |>
  select(
    sample_id,
    treatment_code,
    all_of(soil_vars)
  ) |>
  mutate(
    across(
      all_of(soil_vars),
      ~ as.numeric(scale(.x))
    )
  )

soil_z

write_tsv(
  soil_z,
  result_path("tables", "soil_variables_zscores.tsv")
)

# ------------------------------------------------------------
# 3. Soil-variable heatmap
# ------------------------------------------------------------

# ------------------------------------------------------------
# 4. Pearson correlation matrix among soil variables
# ------------------------------------------------------------

soil_numeric <- env |>
  select(all_of(soil_vars))

soil_cor <- cor(
  soil_numeric,
  method = "pearson",
  use = "pairwise.complete.obs"
)

round(soil_cor, 2)

# ------------------------------------------------------------
# Publication-ready standardized soil heatmap
# ------------------------------------------------------------

soil_heatmap_data <- soil_z |>
  pivot_longer(
    cols = all_of(soil_vars),
    names_to = "variable",
    values_to = "z_score"
  ) |>
  mutate(
    treatment_code = factor(
      treatment_code,
      levels = c("CK", "LS", "MS")
    ),

    sample_id = factor(
      sample_id,
      levels = rev(c(
        "CK-1",
        "CK-2",
        "CK-3",
        "LS-1",
        "LS-2",
        "LS-3",
        "MS-1",
        "MS-2",
        "MS-3"
      ))
    ),

    variable = factor(
      variable,
      levels = c(
        "ec_us_cm",
        "ph",
        "total_n_g_kg",
        "nitrate_n_mg_kg",
        "ammonium_n_mg_kg"
      ),
      labels = c(
        "EC",
        "pH",
        "Total N",
        "NO3-N",
        "NH4-N"
      )
    )
  )

# Symmetric colour limits around zero
z_limit <- max(abs(soil_heatmap_data$z_score))

p_soil_heatmap_pub <- ggplot(
  soil_heatmap_data,
  aes(
    x = variable,
    y = sample_id,
    fill = z_score
  )
) +
  geom_tile(
    colour = "white",
    linewidth = 0.8
  ) +

  geom_text(
    aes(label = sprintf("%.2f", z_score)),
    size = 3.4
  ) +

  facet_grid(
    treatment_code ~ .,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +

  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-z_limit, z_limit),
    name = "Standardized\nvalue (z-score)"
  ) +

  labs(
    x = NULL,
    y = NULL
  ) +

  theme_minimal(base_size = 12) +

  theme(
    panel.grid = element_blank(),

    axis.text.x = element_text(
      colour = "black",
      size = 10
    ),

    axis.text.y = element_text(
      colour = "black",
      size = 10
    ),

    axis.ticks = element_blank(),

    strip.placement = "outside",

    strip.background = element_rect(
      fill = "grey95",
      colour = "grey50",
      linewidth = 0.4
    ),

    strip.text.y.left = element_text(
      angle = 0,
      face = "bold",
      size = 10
    ),

    panel.spacing.y = unit(0.15, "lines"),

    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),

    plot.margin = margin(8, 8, 8, 8)
  )

p_soil_heatmap_pub


ggsave(
  result_path("figures", "standardized_soil_variables_heatmap.png"),
  p_soil_heatmap_pub,
  width = 7.5,
  height = 6,
  units = "in",
  dpi = 600
)

# ------------------------------------------------------------
# Import PCoA coordinates generated in 04_beta_diversity.R
# ------------------------------------------------------------

pcoa_scores <- read_tsv(
  result_path("tables", "pcoa_scores.tsv"),
  show_col_types = FALSE
)

pcoa_variance <- read_tsv(
  result_path("tables", "pcoa_variance_explained.tsv"),
  show_col_types = FALSE
)

pcoa_scores
pcoa_variance


soil_envfit <- env |>
  select(
    sample_id,
    all_of(soil_vars)
  ) |>
  arrange(
    match(sample_id, pcoa_scores$sample)
  )

soil_envfit_matrix <- soil_envfit |>
  select(-sample_id) |>
  as.data.frame()

rownames(soil_envfit_matrix) <- soil_envfit$sample_id

ordination_matrix <- pcoa_scores |>
  select(PCoA1, PCoA2) |>
  as.matrix()

rownames(ordination_matrix) <- pcoa_scores$sample

set.seed(1998)

soil_fit <- envfit(
  ordination_matrix,
  soil_envfit_matrix,
  permutations = 9999
)

soil_fit


# ------------------------------------------------------------
# 5. Extract and save envfit statistics
# ------------------------------------------------------------

soil_envfit_results <- tibble(
  variable = rownames(soil_fit$vectors$arrows),

  PCoA1_direction = soil_fit$vectors$arrows[, 1],
  PCoA2_direction = soil_fit$vectors$arrows[, 2],

  R2 = soil_fit$vectors$r,
  p_value = soil_fit$vectors$pvals
) |>
  mutate(
    p_BH = p.adjust(
      p_value,
      method = "BH"
    ),

    variable = recode(
      variable,
      ec_us_cm = "EC",
      ph = "pH",
      total_n_g_kg = "Total N",
      nitrate_n_mg_kg = "NO3-N",
      ammonium_n_mg_kg = "NH4-N"
    )
  )

soil_envfit_results

write_tsv(
  soil_envfit_results,
  result_path("tables", "soil_envfit_results.tsv")
)
