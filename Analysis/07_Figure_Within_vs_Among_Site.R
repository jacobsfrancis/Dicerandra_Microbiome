# 7) Plot within-site and among-site dissimilarities ####

library(tidyverse)
library(ggdist)

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

# 2) Prepare and draw within-versus-among distributions ####

library(ggdist)

## 2a) Prepare staggered half-violin data ####

within_among_plot_data <- bray_long %>%
  filter(
    is.finite(Bray_Curtis)
  ) %>%
  mutate(
    Compartment_Label = unname(
      compartment_labels[
        Compartment
      ]
    ),

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

    Site_Comparison = factor(
      Same_Site,
      levels = c(
        FALSE,
        TRUE
      ),
      labels = c(
        "Among sites",
        "Within site"
      )
    ),

    Compartment_Number = as.numeric(
      Compartment_Label
    ),

    Plot_Y = Compartment_Number +
      case_when(
        Site_Comparison == "Among sites" ~ 0.18,
        Site_Comparison == "Within site" ~ -0.18
      )
  )

## 2b) Draw and save staggered half-violins ####

within_among_plot <- ggplot(
  within_among_plot_data,
  aes(
    x = Bray_Curtis,
    y = Plot_Y,
    fill = Site_Comparison,
    color = Site_Comparison
  )
) +
  ggdist::stat_halfeye(
    side = "top",
    orientation = "horizontal",
    adjust = 0.7,
    width = 0.34,
    slab_alpha = 0.72,
    slab_linewidth = 0.65,
    point_interval = ggdist::median_qi,
    .width = c(
      0.50,
      0.95
    ),
    point_size = 2.8,
    interval_size = 0.9,
    position = "identity",
    normalize = "groups"
  ) +
  scale_fill_manual(
    values = c(
      "Among sites" = "#D55E00",
      "Within site" = "#0072B2"
    ),
    name = NULL
  ) +
  scale_color_manual(
    values = c(
      "Among sites" = "#D55E00",
      "Within site" = "#0072B2"
    ),
    name = NULL
  ) +
  scale_y_continuous(
    breaks = seq_along(
      levels(
        within_among_plot_data$Compartment_Label
      )
    ),
    labels = levels(
      within_among_plot_data$Compartment_Label
    ),
    expand = expansion(
      mult = c(
        0.04,
        0.08
      )
    )
  ) +
  scale_x_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      by = 0.25
    ),
    expand = expansion(
      mult = c(
        0.01,
        0.02
      )
    )
  ) +
  labs(
    title = paste(
      "Within-site communities are generally more similar",
      "than among-site communities"
    ),
    subtitle = paste(
      "Half-violins show Bray-Curtis dissimilarities;",
      "points are medians and thick/thin lines show",
      "the central 50% and 95% intervals"
    ),
    x = "Bray-Curtis dissimilarity",
    y = NULL,
    caption = paste(
      "All compartments are restricted to the same complete set of sites:",
      paste(
        complete_sites,
        collapse = ", "
      )
    )
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(
        alpha = 0.8
      )
    ),
    color = "none"
  ) +
  theme_classic(
    base_size = 13
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 16
    ),
    plot.subtitle = element_text(
      size = 12
    ),
    axis.title.x = element_text(
      face = "bold",
      size = 13
    ),
    axis.text.y = element_text(
      face = "bold",
      size = 11,
      margin = margin(
        r = 8
      )
    ),
    axis.text.x = element_text(
      size = 11
    ),
    axis.ticks.y = element_blank(),
    legend.position = "top",
    legend.justification = "center",
    legend.key.width = grid::unit(
      1.2,
      "cm"
    ),
    panel.grid.major.x = element_line(
      color = "grey90",
      linewidth = 0.4
    )
  )

ggsave(
  filename = file.path(
    output_directory,
    "within_vs_among_site_half_violins.png"
  ),
  plot = within_among_plot,
  width = 10,
  height = 7,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_directory,
    "within_vs_among_site_half_violins.pdf"
  ),
  plot = within_among_plot,
  width = 10,
  height = 7,
  units = "in"
)
