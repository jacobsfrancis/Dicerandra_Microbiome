# 6) Plot site PERMANOVA effect sizes ####

library(tidyverse)
library(ggtext)

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
nonsignificant_color <- "grey50"

compartment_display_labels <- c(
  Soil_1m  = "Bulk soil (&gt;1 m)",
  Soil_p   = "Rhizosphere soil",
  Root     = "Root epiphyte",
  Leaf_EN  = "Leaf endophyte",
  Leaf_EP  = "Leaf epiphyte",
  Fruit_EN = "Fruit endophyte",
  Fruit_EP = "Fruit epiphyte"
)

color_span <- function(text, color) {
  paste0(
    "<span style='color:", color, ";'>",
    text,
    "</span>"
  )
}

# 3) Prepare site-effect data ####

site_r2_plot_data <- site_permanova_results %>%
  filter(
    Analysis_Status == "Success",
    !is.na(Site_R2),
    is.finite(Site_R2)
  )

if (nrow(site_r2_plot_data) == 0) {
  stop(
    paste(
      "No successful site PERMANOVA results were available to plot.",
      "Inspect site_permanova_results$Analysis_Status."
    )
  )
}

# Recover compartment codes from readable labels if the results table does not
# already contain a Compartment column.
if (!"Compartment" %in% names(site_r2_plot_data)) {
  site_r2_plot_data <- site_r2_plot_data %>%
    mutate(
      Compartment = case_when(
        Compartment_Label == "Bulk soil (>1 m)"  ~ "Soil_1m",
        Compartment_Label == "Rhizosphere soil" ~ "Soil_p",
        Compartment_Label == "Root epiphyte"     ~ "Root",
        Compartment_Label == "Root surface"      ~ "Root",
        Compartment_Label == "Leaf endophyte"    ~ "Leaf_EN",
        Compartment_Label == "Leaf epiphyte"     ~ "Leaf_EP",
        Compartment_Label == "Fruit endophyte"   ~ "Fruit_EN",
        Compartment_Label == "Fruit epiphyte"    ~ "Fruit_EP",
        TRUE                                       ~ NA_character_
      )
    )
}

if (any(is.na(site_r2_plot_data$Compartment))) {
  stop(
    paste(
      "At least one plotted row could not be matched to a compartment color.",
      "Inspect site_r2_plot_data$Compartment and Compartment_Label."
    )
  )
}

# Preserve the original y-axis order while replacing labels with colored HTML.
if (is.factor(site_r2_plot_data$Compartment_Label)) {
  original_label_order <- levels(site_r2_plot_data$Compartment_Label)
} else {
  original_label_order <- unique(site_r2_plot_data$Compartment_Label)
}

label_lookup <- site_r2_plot_data %>%
  distinct(Compartment_Label, Compartment) %>%
  mutate(
    Colored_Label = map2_chr(
      compartment_display_labels[Compartment],
      compartment_colors[Compartment],
      color_span
    )
  )

colored_label_order <- label_lookup$Colored_Label[
  match(original_label_order, label_lookup$Compartment_Label)
]
colored_label_order <- colored_label_order[!is.na(colored_label_order)]

# Add the R-squared label just to the right of each point, then calculate an
# x-axis limit that leaves enough room for the complete label.
observed_range <- max(site_r2_plot_data$Site_R2, na.rm = TRUE)
label_offset <- max(0.012, observed_range * 0.045)
right_padding <- max(0.065, observed_range * 0.20)

site_r2_plot_data <- site_r2_plot_data %>%
  left_join(label_lookup, by = c("Compartment_Label", "Compartment")) %>%
  mutate(
    Colored_Label = factor(
      Colored_Label,
      levels = colored_label_order
    ),
    R2_Label = sprintf("R² = %.2f", Site_R2),
    R2_Label_X = Site_R2 + label_offset
  )

site_r2_upper_limit <- max(
  0.10,
  max(site_r2_plot_data$R2_Label_X, na.rm = TRUE) + right_padding
)

# 4) Plot ####

site_r2_plot <- ggplot(
  site_r2_plot_data,
  aes(
    x = Site_R2,
    y = Colored_Label,
    color = PERMANOVA_Result,
    shape = Dispersion_Result
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = Site_R2,
      yend = Colored_Label
    ),
    color = "grey75",
    linewidth = 0.8
  ) +
  geom_point(
    size = 4.5
  ) +
  geom_text(
    aes(
      label = sprintf("R² = %.2f", Site_R2)
    ),
    nudge_x = -0.0,   # slightly left of point
    nudge_y = 0.18,    # above line
    hjust = 1,
    vjust = 0,
    color = "#995C8F",
    fontface = "bold",
    size = 7
  ) +
  scale_color_manual(
    values = c(
      "Site effect: BH-adjusted p < 0.05" = significant_color,
      "Site effect: BH-adjusted p >= 0.05" = nonsignificant_color
    ),
    name = NULL
  ) +
  scale_shape_manual(
    values = c(
      "No significant dispersion difference" = 16,
      "Dispersion differs among sites" = 17
    ),
    name = "Multivariate dispersion"
  ) +
  scale_x_continuous(
    limits = c(0, site_r2_upper_limit),
    expand = expansion(mult = c(0, 0))
  ) +
  coord_cartesian(
    clip = "off"
  ) +
  labs(
    x = expression("Site PERMANOVA " * R^2),
    y = NULL
  ) +
  theme_classic(
    base_size = 20,
    base_family = "Helvetica"
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),

    axis.title.x = element_text(
      size = 20,
      face = "bold"
    ),

    axis.text.y = ggtext::element_markdown(
      size = 20,
      lineheight = 1.05,
      margin = margin(r = 8),
      face="bold"
    ),

    axis.ticks.y = element_blank(),
    legend.position = "none",
    panel.background = element_blank(),
    plot.background = element_blank(),

    plot.margin = margin(
      t = 6,
      r = 18,
      b = 6,
      l = 6
    )
  )

# 5) Save ####

figure_width <- 6.5
figure_height <- 9.5

ggsave(
  filename = file.path(
    output_directory,
    "site_PERMANOVA_R2_by_compartment.png"
  ),
  plot = site_r2_plot,
  width = figure_width,
  height = figure_height,
  units = "in",
  dpi = 600,
  bg = "transparent"
)

ggsave(
  filename = file.path(
    output_directory,
    "site_PERMANOVA_R2_by_compartment.pdf"
  ),
  plot = site_r2_plot,
  width = figure_width,
  height = figure_height,
  units = "in",
  bg = "transparent"
)
