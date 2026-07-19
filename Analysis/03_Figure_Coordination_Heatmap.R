# 3) Plot the compartment-coordination heat map ####

library(tidyverse)

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

# 2) Prepare and draw the coordination heat map ####

## 2a) Prepare lower-triangle heat-map data ####

heatmap_data <- distance_model_results %>%
  filter(
    Model_Status == "Success"
  ) %>%
  mutate(
    Predictor = factor(
      Predictor,
      levels = compartment_order
    ),

    Response = factor(
      Response,
      levels = rev(
        compartment_order
      )
    ),

    Estimate_Label = if_else(
      is.na(Estimate),
      "",
      sprintf(
        "%.2f",
        Estimate
      )
    ),

    Significance_Label = case_when(
      !is.na(P_Adjusted_BH) &
        P_Adjusted_BH < 0.001 ~
        "***",

      !is.na(P_Adjusted_BH) &
        P_Adjusted_BH < 0.01 ~
        "**",

      !is.na(P_Adjusted_BH) &
        P_Adjusted_BH < 0.05 ~
        "*",

      TRUE ~ ""
    ),

    Tile_Label = paste0(
      Estimate_Label,
      Significance_Label
    )
  )

heatmap_limit <- max(
  abs(
    heatmap_data$Estimate
  ),
  na.rm = TRUE
)

if (
  !is.finite(heatmap_limit) ||
  heatmap_limit == 0
) {
  heatmap_limit <- 1
}

coefficient_heatmap <- ggplot(
  heatmap_data,
  aes(
    x = Predictor,
    y = Response,
    fill = Estimate
  )
) +
  geom_tile(
    color = "white",
    linewidth = 0.8
  ) +
  geom_text(
    aes(
      label = Tile_Label
    ),
    size = 3.8,
    fontface = "bold"
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(
      -heatmap_limit,
      heatmap_limit
    ),
    name =
      "Standardized\ncoefficient"
  ) +
  scale_x_discrete(
    labels =
      compartment_labels,
    drop = FALSE
  ) +
  scale_y_discrete(
    labels =
      compartment_labels,
    drop = FALSE
  ) +
  coord_fixed() +
  labs(
    title =
      paste(
        "Coordination of microbiome",
        "dissimilarity among compartments"
      ),

    subtitle =
      paste(
        "All plant pairs;",
        "models control Same_Site",
        "and include crossed site",
        "and plant random intercepts"
      ),

    x = "Predictor compartment",
    y = "Response compartment",

    caption =
      paste(
        "* q < 0.05;",
        "** q < 0.01;",
        "*** q < 0.001"
      )
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    panel.grid =
      element_blank(),

    axis.text.x =
      element_text(
        angle = 45,
        hjust = 1,
        vjust = 1
      ),

    axis.title =
      element_text(
        face = "bold"
      ),

    plot.title =
      element_text(
        face = "bold"
      ),

    legend.title =
      element_text(
        face = "bold"
      )
  )

ggsave(
  filename = file.path(
    output_directory,
    "all_pair_coefficient_heatmap.png"
  ),
  plot = coefficient_heatmap,
  width = 10,
  height = 8,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_directory,
    "all_pair_coefficient_heatmap.pdf"
  ),
  plot = coefficient_heatmap,
  width = 10,
  height = 8,
  units = "in"
)
