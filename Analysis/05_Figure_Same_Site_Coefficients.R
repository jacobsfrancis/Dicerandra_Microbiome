# 5) Plot same-site coefficients from coordination models ####

library(tidyverse)
library(forcats)

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

# 2) Prepare and draw same-site coefficient plot ####

## 2a) Prepare same-site coefficient data ####

same_site_plot_data <-
  distance_model_results %>%
  filter(
    Model_Status == "Success",
    is.finite(Same_Site_Estimate),
    is.finite(Same_Site_CI_Lower),
    is.finite(Same_Site_CI_Upper)
  ) %>%
  mutate(
    Predictor_Label =
      unname(
        compartment_labels[
          Predictor
        ]
      ),

    Response_Label =
      unname(
        compartment_labels[
          Response
        ]
      ),

    Comparison_Label =
      paste0(
        Response_Label,
        "  ~  ",
        Predictor_Label
      ),

    FDR_Result = if_else(
      Same_Site_Significant_BH,
      "BH-adjusted p < 0.05",
      "BH-adjusted p >= 0.05"
    ),

    Comparison_Category =
      factor(
        Comparison_Category,
        levels =
          comparison_category_order
      )
  ) %>%
  arrange(
    Comparison_Category,
    Same_Site_Estimate
  ) %>%
  mutate(
    Comparison_Label =
      forcats::fct_inorder(
        Comparison_Label
      )
  )

same_site_plot <- ggplot(
  same_site_plot_data,
  aes(
    x = Same_Site_Estimate,
    y = Comparison_Label,
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
      xmin =
        Same_Site_CI_Lower,
      xmax =
        Same_Site_CI_Upper
    ),
    height = 0.18,
    linewidth = 0.8
  ) +
  geom_point(
    size = 2.8
  ) +
  facet_grid(
    rows =
      vars(
        Comparison_Category
      ),
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +
  scale_color_manual(
    values = c(
      "BH-adjusted p < 0.05" =
        "#B2182B",

      "BH-adjusted p >= 0.05" =
        "#4D4D4D"
    ),
    name = NULL
  ) +
  labs(
    title =
      paste(
        "Same-site effect in",
        "all-pair dissimilarity models"
      ),

    subtitle =
      paste(
        "Negative estimates mean",
        "within-site pairs have lower",
        "response dissimilarity than",
        "between-site pairs"
      ),

    x =
      "Same-site coefficient (95% CI)",

    y = NULL
  ) +
  theme_minimal(
    base_size = 12
  ) +
  theme(
    panel.grid.major.y =
      element_blank(),

    panel.grid.minor =
      element_blank(),

    axis.title.x =
      element_text(
        face = "bold"
      ),

    axis.text.y =
      element_text(
        size = 10
      ),

    plot.title =
      element_text(
        face = "bold"
      ),

    legend.position =
      "top",

    strip.placement =
      "outside",

    strip.text.y.left =
      element_text(
        angle = 0,
        face = "bold",
        hjust = 1
      ),

    strip.background =
      element_rect(
        fill = "grey92",
        color = NA
      ),

    panel.spacing.y =
      grid::unit(
        0.7,
        "lines"
      )
  )

ggsave(
  filename = file.path(
    output_directory,
    "same_site_effect_plot_95CI.png"
  ),
  plot = same_site_plot,
  width = 12,
  height = 11,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_directory,
    "same_site_effect_plot_95CI.pdf"
  ),
  plot = same_site_plot,
  width = 12,
  height = 11,
  units = "in"
)

## 3) Summarize model diagnostics ####

model_diagnostics <-
  distance_model_results %>%
  summarise(
    N_Models = n(),

    N_Successful =
      sum(
        Model_Status == "Success"
      ),

    N_Singular =
      sum(
        Singular %in% TRUE,
        na.rm = TRUE
      ),

    N_With_Convergence_Message =
      sum(
        !is.na(
          Convergence_Message
        ) &
          Convergence_Message != ""
      ),

    N_Raw_Significant =
      sum(
        Significant_Raw,
        na.rm = TRUE
      ),

    N_BH_Significant =
      sum(
        Significant_BH,
        na.rm = TRUE
      ),

    N_Same_Site_Raw_Significant =
      sum(
        Same_Site_Significant_Raw,
        na.rm = TRUE
      ),

    N_Same_Site_BH_Significant =
      sum(
        Same_Site_Significant_BH,
        na.rm = TRUE
      )
  )

write_csv(
  model_diagnostics,
  file.path(
    output_directory,
    "model_diagnostic_summary.csv"
  )
)
