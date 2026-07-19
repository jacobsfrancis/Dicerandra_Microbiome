# 10) Plot among-site minus within-site dissimilarity gaps ####

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

# 2) Prepare and draw the compartment-specific site-gap plot ####

site_gap_plot_data <- site_gap_results %>%
  mutate(
    Compartment_Label = factor(
      Compartment_Label,
      levels = rev(
        unname(
          compartment_labels[
            compartment_order
          ]
        )
      )
    ),

    FDR_Result = if_else(
      Significant_BH,
      "BH-adjusted p < 0.05",
      "BH-adjusted p >= 0.05"
    )
  )

site_gap_plot <- ggplot(
  site_gap_plot_data,
  aes(
    x = Among_Minus_Within,
    y = Compartment_Label,
    color = FDR_Result
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.7,
    color = "grey40"
  ) +
  geom_errorbarh(
    aes(
      xmin = CI_Lower,
      xmax = CI_Upper
    ),
    height = 0.16,
    linewidth = 0.9
  ) +
  geom_point(
    size = 3.5
  ) +
  scale_color_manual(
    values = c(
      "BH-adjusted p < 0.05" = "#B2182B",
      "BH-adjusted p >= 0.05" = "#4D4D4D"
    ),
    name = NULL
  ) +
  labs(
    title = "Site-associated similarity differs among microbiome compartments",

    subtitle = paste(
      "Points show the estimated increase in Bray-Curtis dissimilarity",
      "for among-site relative to within-site plant pairs"
    ),

    x = "Among-site minus within-site dissimilarity (95% CI)",
    y = NULL,

    caption = paste(
      "Larger positive estimates indicate stronger site structuring;",
      "models account for repeated plant and plant-pair observations."
    )
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),

    axis.title.x = element_text(
      face = "bold"
    ),

    axis.text.y = element_text(
      face = "bold"
    ),

    legend.position = "top"
  )

ggsave(
  filename = file.path(
    output_directory,
    "among_minus_within_site_gap_by_compartment.png"
  ),
  plot = site_gap_plot,
  width = 9,
  height = 6,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_directory,
    "among_minus_within_site_gap_by_compartment.pdf"
  ),
  plot = site_gap_plot,
  width = 9,
  height = 6,
  units = "in"
)
