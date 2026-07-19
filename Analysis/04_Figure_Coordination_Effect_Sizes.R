# 4) Plot compartment-coordination effect sizes ####

library(tidyverse)
library(forcats)
library(ggtext)
library(ggh4x)

# 0) File paths ####

data_directory <- "../Data/"
output_directory <- file.path("../Output/", "Poster_Common_5_Sites")
intermediate_directory <- file.path(output_directory, "Intermediate")

dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
dir.create(intermediate_directory, showWarnings = FALSE, recursive = TRUE)

# 1) Load shared data and analysis results ####

clean_data_objects <- readRDS(
  file.path(intermediate_directory, "clean_data_objects.rds")
)

analysis_results <- readRDS(
  file.path(intermediate_directory, "analysis_results.rds")
)

list2env(clean_data_objects, envir = .GlobalEnv)
list2env(analysis_results, envir = .GlobalEnv)

# 2) Visual vocabulary ####

compartment_colors <- c(
  Soil_1m  = "#C83A06", # bulk soil
  Soil_p   = "#DD6038", # rhizosphere
  Root     = "#F2C54D", # root surface
  Leaf_EN  = "#476C35", # leaf endophyte
  Leaf_EP  = "#B6CD78", # leaf epiphyte
  Fruit_EN = "#355EAA", # fruit endophyte
  Fruit_EP = "#6EA4CE"  # fruit epiphyte
)

significant_color <- "#995C8F"
nonsignificant_color <- "grey55"

# Pale pathway backgrounds and matching strip-label colors.
category_fill_colors <- c(
  "Belowground pathway" = "#F9E8DF",
  "Leaf pathway"        = "#EDF3DF",
  "Fruit pathway"       = "#E5EEF7",
  "Belowground-leaf"    = "#F1EAF0",
  "Belowground-fruit"   = "#F1EAF0",
  "Belowground−fruit"   = "#F1EAF0",
  "Leaf−fruit"          = "#F1EAF0",
  "Leaf-fruit"          = "#F1EAF0"
)

category_text_colors <- c(
  "Belowground pathway" = "#C83A06",
  "Leaf pathway"        = "#476C35",
  "Fruit pathway"       = "#355EAA",
  "Belowground−leaf"    = "#995C8F",
  "Belowground-leaf"    = "#995C8F",
  "Belowground−fruit"   = "#995C8F",
  "Belowground-fruit"   = "#995C8F",
  "Leaf−fruit"          = "#995C8F",
  "Leaf-fruit"          = "#995C8F"
)

theme_set(
  theme_classic(
    base_size = 18,
    base_family = "Helvetica"
  ) +
    theme(
      axis.title = element_text(
      ),
      plot.title = element_text(
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 2
        )
      ),
      legend.title = element_blank(),
      legend.key = element_blank(),
      panel.background = element_blank(),
      plot.background = element_blank(),
      legend.background = element_blank(),
      plot.margin = margin(
        t = 0,
        r = 0,
        b = 0,
        l = 0
      )
    )
)


# 3) Helpers for colored labels ####

color_span <- function(text, color) {
  paste0(
    "<span style='color:", color, ";'>",
    text,
    "</span>"
  )
}

make_comparison_label <- function(response, predictor) {
  response_text <- unname(compartment_labels[response])
  predictor_text <- unname(compartment_labels[predictor])

  paste0(
    color_span(response_text, compartment_colors[response]),
    " <span style='color:#4D4D4D;'>~</span> ",
    color_span(predictor_text, compartment_colors[predictor])
  )
}

# 4) Prepare categorized coefficient data ####

coefficient_plot_data <-
  distance_model_results %>%
  filter(
    Model_Status == "Success",
    is.finite(Estimate),
    is.finite(CI_Lower),
    is.finite(CI_Upper)
  ) %>%
  mutate(
    Comparison_Label = map2_chr(
      Response,
      Predictor,
      make_comparison_label
    ),

    Significance = if_else(
      Significant_BH,
      "Significant",
      "Not significant"
    ),

    Comparison_Category = factor(
      Comparison_Category,
      levels = comparison_category_order
    )
  ) %>%
  arrange(
    Comparison_Category,
    Estimate
  ) %>%
  mutate(
    Comparison_Label = fct_inorder(Comparison_Label)
  )

# Keep strip styling synchronized with the category order actually present.
category_levels <- levels(droplevels(coefficient_plot_data$Comparison_Category))

strip_fill <- unname(category_fill_colors[category_levels])
strip_text <- unname(category_text_colors[category_levels])

# Fallbacks prevent errors if a category label differs slightly.
strip_fill[is.na(strip_fill)] <- "grey95"
strip_text[is.na(strip_text)] <- "grey25"

# 5) Plot ####

coefficient_plot <- ggplot(
  coefficient_plot_data,
  aes(
    x = Estimate,
    y = Comparison_Label,
    color = Significance,
    size = Significance
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.65,
    color = "grey45"
  ) +
  geom_errorbarh(
    aes(
      xmin = CI_Lower,
      xmax = CI_Upper
    ),
    height = 0.16,
    linewidth = 0.8
  ) +
  geom_point() +
  ggh4x::facet_grid2(
    rows = vars(Comparison_Category),
    scales = "free_y",
    space = "free_y",
    switch = "y",
    axes = "all",
    remove_labels = "x",
    strip = ggh4x::strip_themed(
      background_y = lapply(
        strip_fill,
        function(x) element_rect(fill = x, color = NA)
      ),
      text_y = Map(
        function(color_value) {
          element_text(
            color = color_value,
            face = "bold",
            angle = 0,
            hjust = 1,
            size=22
          )
        },
        strip_text
      )
    )
  ) +
  scale_color_manual(
    values = c(
      "Significant" = significant_color,
      "Not significant" = nonsignificant_color
    )
  ) +
  scale_size_manual(
    values = c(
      "Significant" = 3.5,
      "Not significant" = 2.4
    )
  ) +
  labs(
    x = "Standardized coefficient (95% CI)",
    y = NULL
  ) +
  theme_classic(
    base_size = 18,
    base_family = "Helvetica"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),

    axis.title.x = element_text(
      face = "bold",
      margin = margin(t = 8),
      size=20
    ),

    axis.text.y = ggtext::element_markdown(
      size = 20,
      color = "grey20",
      lineheight = 1.05
    ),

    legend.position = "none",

    strip.placement = "outside",

    panel.spacing.y = grid::unit(0.12, "lines"),

    panel.background = element_blank(),
    plot.background = element_blank(),

    plot.margin = margin(
      t = 4,
      r = 8,
      b = 4,
      l = 4
    )
  )

# 6) Save ####

ggsave(
  filename = file.path(
    output_directory,
    "all_pair_coefficient_plot_95CI_by_category.png"
  ),
  plot = coefficient_plot,
  width = 14,
  height = 9.5,
  units = "in",
  dpi = 600,
  bg = "transparent"
)

ggsave(
  filename = file.path(
    output_directory,
    "all_pair_coefficient_plot_95CI_by_category.pdf"
  ),
  plot = coefficient_plot,
  width = 14,
  height = 9.5,
  units = "in",
  bg = "transparent"
)
