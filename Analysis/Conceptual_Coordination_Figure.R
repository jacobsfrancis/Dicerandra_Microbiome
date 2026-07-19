# 10) Conceptual figure: microbiome coordination ####
#
# Run this section after the data-cleaning portion of the main analysis has
# created:
#
#   metadata_matched
#   taxa_abundance_matched
#   bray_pair_table
#   output_directory
#
# Panel A shows bulk-soil microbiomes (Soil_1m).
# Panel B shows rhizosphere microbiomes (Soil_p).
# Panel C shows the Bray-Curtis dissimilarities of all available plant pairs.
#
# Colors identify microbiome compartments, not sampling sites.
# All observations in Panel C are charcoal.

library(tidyverse)
library(vegan)
library(patchwork)
library(ggrepel)

## 10a) Define the compartment and neutral colors ####

compartment_colors <- c(
  Soil_1m  = "#C83A06", # bulk soil: red
  Soil_p   = "#DD6038", # rhizosphere: orange
  Root     = "#F2C54D", # root: yellow
  Leaf_EN  = "#476C35", # leaf endophyte: dark green
  Leaf_EP  = "#B6CD78", # leaf epiphyte: light green
  Fruit_EN = "#355EAA", # fruit endophyte: dark blue
  Fruit_EP = "#6EA4CE"  # fruit epiphyte: light blue
)

figure_colors <- c(
  charcoal = "#262626",
  gray = "#777777",
  light_gray = "#D8D8D8"
)

## 10b) Calculate PCoA ordinations ####

make_compartment_pcoa <- function(
    compartment_name,
    metadata_table,
    abundance_matrix
) {

  compartment_metadata <- metadata_table %>%
    filter(
      Compartment == compartment_name
    ) %>%
    arrange(
      Site,
      suppressWarnings(as.numeric(Plant)),
      Plant_ID
    )

  compartment_abundance <- abundance_matrix[
    compartment_metadata$Sample_Index,
    ,
    drop = FALSE
  ]

  keep_samples <- rowSums(
    compartment_abundance,
    na.rm = TRUE
  ) > 0

  compartment_metadata <- compartment_metadata[
    keep_samples,
    ,
    drop = FALSE
  ]

  compartment_abundance <- compartment_abundance[
    keep_samples,
    ,
    drop = FALSE
  ]

  keep_taxa <- colSums(
    compartment_abundance,
    na.rm = TRUE
  ) > 0

  compartment_abundance <- compartment_abundance[
    ,
    keep_taxa,
    drop = FALSE
  ]

  if (nrow(compartment_abundance) < 3) {
    stop(
      paste(
        "Fewer than three nonempty samples were available for",
        compartment_name
      )
    )
  }

  if (ncol(compartment_abundance) < 1) {
    stop(
      paste(
        "No taxa remained for",
        compartment_name
      )
    )
  }

  bray_distance <- vegan::vegdist(
    compartment_abundance,
    method = "bray"
  )

  pcoa_fit <- stats::cmdscale(
    bray_distance,
    k = 2,
    eig = TRUE,
    add = TRUE
  )

  score_table <- as.data.frame(
    pcoa_fit$points
  )

  names(score_table) <- c(
    "Axis_1",
    "Axis_2"
  )

  score_table$Sample_Index <- rownames(
    score_table
  )

  score_table <- score_table %>%
    as_tibble() %>%
    left_join(
      compartment_metadata %>%
        transmute(
          Sample_Index = as.character(Sample_Index),
          Site,
          Plant,
          Plant_ID,
          Compartment
        ),
      by = "Sample_Index"
    )

  positive_eigenvalues <- pcoa_fit$eig[
    pcoa_fit$eig > 0
  ]

  variance_explained <- 100 *
    positive_eigenvalues /
    sum(positive_eigenvalues)

  list(
    scores = score_table,
    axis_1_percent = variance_explained[1],
    axis_2_percent = variance_explained[2]
  )
}

bulk_pcoa <- make_compartment_pcoa(
  compartment_name = "Soil_1m",
  metadata_table = metadata_matched,
  abundance_matrix = taxa_abundance_matched
)

rhizosphere_pcoa <- make_compartment_pcoa(
  compartment_name = "Soil_p",
  metadata_table = metadata_matched,
  abundance_matrix = taxa_abundance_matched
)

## 10c) Select one representative within-site pair ####
#
# The pair is selected from plants represented in both ordinations.
# It is chosen near the 70th percentile of Bray-Curtis dissimilarity in
# both compartments so it is visible without being an extreme outlier.
#
# To choose a pair manually, replace the automatic-selection block with:
#
# example_plant_1 <- "2_1"
# example_plant_2 <- "2_2"
#
# example_pair <- bray_pair_table %>%
#   filter(
#     (Plant_1_ID == example_plant_1 &
#        Plant_2_ID == example_plant_2) |
#     (Plant_1_ID == example_plant_2 &
#        Plant_2_ID == example_plant_1)
#   ) %>%
#   slice(1)

plants_in_both_ordinations <- intersect(
  bulk_pcoa$scores$Plant_ID,
  rhizosphere_pcoa$scores$Plant_ID
)

example_pair_candidates <- bray_pair_table %>%
  filter(
    Same_Site,
    Plant_1_ID %in% plants_in_both_ordinations,
    Plant_2_ID %in% plants_in_both_ordinations,
    is.finite(Soil_1m),
    is.finite(Soil_p)
  ) %>%
  mutate(
    Bulk_Percentile = percent_rank(Soil_1m),
    Rhizosphere_Percentile = percent_rank(Soil_p),
    Selection_Score =
      abs(Bulk_Percentile - 0.70) +
      abs(Rhizosphere_Percentile - 0.70)
  ) %>%
  arrange(
    Selection_Score
  )

if (nrow(example_pair_candidates) == 0) {
  stop(
    paste(
      "No within-site pair was represented in both the bulk-soil",
      "and rhizosphere ordinations."
    )
  )
}

example_pair <- example_pair_candidates %>%
  slice(1)

example_plant_1 <- example_pair$Plant_1_ID
example_plant_2 <- example_pair$Plant_2_ID
example_site <- example_pair$Site_1

## 10d) Prepare the highlighted pair in each ordination ####

prepare_example_pair <- function(
    score_table,
    plant_1,
    plant_2,
    distance_value
) {

  highlighted_points <- score_table %>%
    filter(
      Plant_ID %in% c(
        plant_1,
        plant_2
      )
    ) %>%
    mutate(
      Example_Label = case_when(
        Plant_ID == plant_1 ~ "Plant A",
        Plant_ID == plant_2 ~ "Plant B",
        TRUE ~ NA_character_
      )
    ) %>%
    arrange(
      factor(
        Example_Label,
        levels = c(
          "Plant A",
          "Plant B"
        )
      )
    )

  if (nrow(highlighted_points) != 2) {
    stop(
      paste(
        "The selected pair was not represented exactly once",
        "in one of the ordinations."
      )
    )
  }

  connecting_line <- tibble(
    x = highlighted_points$Axis_1[1],
    y = highlighted_points$Axis_2[1],
    xend = highlighted_points$Axis_1[2],
    yend = highlighted_points$Axis_2[2],
    label_x = mean(highlighted_points$Axis_1),
    label_y = mean(highlighted_points$Axis_2),
    Distance_Label = sprintf(
      "Bray-Curtis = %.2f",
      distance_value
    )
  )

  list(
    points = highlighted_points,
    line = connecting_line
  )
}

bulk_example <- prepare_example_pair(
  score_table = bulk_pcoa$scores,
  plant_1 = example_plant_1,
  plant_2 = example_plant_2,
  distance_value = example_pair$Soil_1m
)

rhizosphere_example <- prepare_example_pair(
  score_table = rhizosphere_pcoa$scores,
  plant_1 = example_plant_1,
  plant_2 = example_plant_2,
  distance_value = example_pair$Soil_p
)

## 10e) Draw the ordination panels ####

make_ordination_panel <- function(
    pcoa_object,
    example_object,
    panel_title,
    compartment_color,
    panel_letter
) {

  ggplot(
    pcoa_object$scores,
    aes(
      x = Axis_1,
      y = Axis_2
    )
  ) +

    geom_point(
      shape = 21,
      size = 3.0,
      stroke = 0.35,
      fill = scales::alpha(
        compartment_color,
        0.30
      ),
      color = scales::alpha(
        compartment_color,
        0.48
      )
    ) +

    geom_segment(
      data = example_object$line,
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ),
      inherit.aes = FALSE,
      color = figure_colors[["charcoal"]],
      linewidth = 1.05,
      linetype = "dashed"
    ) +

    geom_point(
      data = example_object$points,
      aes(
        x = Axis_1,
        y = Axis_2
      ),
      inherit.aes = FALSE,
      shape = 21,
      size = 5.5,
      stroke = 1.05,
      fill = compartment_color,
      color = figure_colors[["charcoal"]]
    ) +

    ggrepel::geom_text_repel(
      data = example_object$points,
      aes(
        x = Axis_1,
        y = Axis_2,
        label = Example_Label
      ),
      inherit.aes = FALSE,
      size = 3.8,
      fontface = "bold",
      color = figure_colors[["charcoal"]],
      box.padding = 0.45,
      point.padding = 0.45,
      min.segment.length = 0,
      segment.color = figure_colors[["gray"]]
    ) +

    ggrepel::geom_label_repel(
      data = example_object$line,
      aes(
        x = label_x,
        y = label_y,
        label = Distance_Label
      ),
      inherit.aes = FALSE,
      size = 3.5,
      fontface = "bold",
      color = compartment_color,
      fill = scales::alpha(
        "white",
        0.95
      ),
      label.size = 0,
      box.padding = 0.35,
      point.padding = 0.25,
      min.segment.length = 0,
      segment.color = figure_colors[["gray"]]
    ) +

    labs(
      title = paste0(
        panel_letter,
        ". ",
        panel_title
      ),
      subtitle = paste0(
        "The same two plants are highlighted from site ",
        example_site
      ),
      x = sprintf(
        "PCoA 1 (%.1f%%)",
        pcoa_object$axis_1_percent
      ),
      y = sprintf(
        "PCoA 2 (%.1f%%)",
        pcoa_object$axis_2_percent
      )
    ) +

    coord_equal() +

    theme_classic(
      base_size = 12
    ) +

    theme(
      plot.title = element_text(
        face = "bold",
        size = 14
      ),
      plot.subtitle = element_text(
        size = 10,
        color = figure_colors[["gray"]]
      ),
      axis.title = element_text(
        face = "bold"
      ),
      axis.text = element_text(
        color = figure_colors[["charcoal"]]
      )
    )
}

bulk_ordination_plot <- make_ordination_panel(
  pcoa_object = bulk_pcoa,
  example_object = bulk_example,
  panel_title = "Bulk-soil microbiomes",
  compartment_color = compartment_colors[["Soil_1m"]],
  panel_letter = "A"
)

rhizosphere_ordination_plot <- make_ordination_panel(
  pcoa_object = rhizosphere_pcoa,
  example_object = rhizosphere_example,
  panel_title = "Rhizosphere microbiomes",
  compartment_color = compartment_colors[["Soil_p"]],
  panel_letter = "B"
)

## 10f) Prepare the all-pair coordination scatterplot ####

coordination_scatter_data <- bray_pair_table %>%
  filter(
    is.finite(Soil_1m),
    is.finite(Soil_p)
  ) %>%
  mutate(
    Pair_Type = factor(
      if_else(
        Same_Site,
        "Within site",
        "Among sites"
      ),
      levels = c(
        "Among sites",
        "Within site"
      )
    ),
    Is_Example_Pair =
      (
        Plant_1_ID == example_plant_1 &
          Plant_2_ID == example_plant_2
      ) |
      (
        Plant_1_ID == example_plant_2 &
          Plant_2_ID == example_plant_1
      )
  )

## 10g) Draw the coordination scatterplot ####
#
# Every pair is charcoal. Alpha distinguishes within-site and among-site
# pairs without assigning colors to sampling sites.
#
# The fitted line is only a visual guide. The inferential mixed model uses
# all pairs and includes Same_Site plus crossed site and plant random effects.

coordination_scatter_plot <- ggplot(
  coordination_scatter_data,
  aes(
    x = Soil_1m,
    y = Soil_p
  )
) +

  geom_point(
    data = coordination_scatter_data %>%
      filter(
        !Same_Site
      ),
    color = figure_colors[["charcoal"]],
    size = 1.8,
    alpha = 0.15
  ) +

  geom_point(
    data = coordination_scatter_data %>%
      filter(
        Same_Site
      ),
    color = figure_colors[["charcoal"]],
    size = 2.1,
    alpha = 0.50
  ) +

  geom_smooth(
    method = "lm",
    formula = y ~ x,
    se = TRUE,
    color = figure_colors[["charcoal"]],
    fill = figure_colors[["light_gray"]],
    linewidth = 1.15,
    alpha = 0.28
  ) +

  geom_point(
    data = coordination_scatter_data %>%
      filter(
        Is_Example_Pair
      ),
    shape = 21,
    size = 5.5,
    stroke = 1.35,
    fill = "white",
    color = figure_colors[["charcoal"]]
  ) +

  geom_point(
    data = coordination_scatter_data %>%
      filter(
        Is_Example_Pair
      ),
    size = 2.3,
    color = figure_colors[["charcoal"]]
  ) +

  ggrepel::geom_label_repel(
    data = coordination_scatter_data %>%
      filter(
        Is_Example_Pair
      ),
    aes(
      label = "Plant A-B pair"
    ),
    size = 3.5,
    fontface = "bold",
    color = figure_colors[["charcoal"]],
    fill = "white",
    label.size = 0.25,
    box.padding = 0.45,
    point.padding = 0.45,
    min.segment.length = 0,
    segment.color = figure_colors[["charcoal"]]
  ) +

  annotate(
    geom = "label",
    x = -Inf,
    y = Inf,
    hjust = -0.05,
    vjust = 1.10,
    label = paste(
      "Each point = one plant pair",
      "Darker = within-site pair",
      "Lighter = among-site pair",
      sep = "\n"
    ),
    size = 3.3,
    color = figure_colors[["charcoal"]],
    fill = scales::alpha(
      "white",
      0.92
    ),
    label.size = 0
  ) +

  labs(
    title = "C. Coordination of microbiome dissimilarity",
    subtitle = paste(
      "Do plant pairs that differ more in bulk soil",
      "also differ more in their rhizosphere?"
    ),
    x = "Bulk-soil Bray-Curtis dissimilarity",
    y = "Rhizosphere Bray-Curtis dissimilarity",
    caption = paste(
      "The line is a visual guide. The inferential mixed model uses all",
      "plant pairs and accounts for same-site status and repeated use",
      "of sites and plants."
    )
  ) +

  theme_classic(
    base_size = 12
  ) +

  theme(
    plot.title = element_text(
      face = "bold",
      size = 14
    ),
    plot.subtitle = element_text(
      size = 10,
      color = figure_colors[["gray"]]
    ),
    axis.title = element_text(
      face = "bold"
    ),
    axis.text = element_text(
      color = figure_colors[["charcoal"]]
    ),
    plot.caption = element_text(
      hjust = 0,
      size = 9,
      color = figure_colors[["gray"]]
    ),
    legend.position = "none"
  )

## 10h) Assemble and save the conceptual figure ####

coordination_method_figure <-
  (
    bulk_ordination_plot |
      rhizosphere_ordination_plot
  ) /
    coordination_scatter_plot +

  patchwork::plot_layout(
    heights = c(
      1,
      1.08
    )
  ) +

  patchwork::plot_annotation(
    title =
      "How we tested coordination of microbiome turnover across compartments",
    subtitle = paste(
      "For every pair of plants, we calculated microbiome dissimilarity",
      "in each compartment and tested whether those distances covaried."
    ),
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 18,
        hjust = 0.5
      ),
      plot.subtitle = element_text(
        size = 12,
        hjust = 0.5
      )
    )
  )

coordination_method_figure

ggsave(
  filename = file.path(
    output_directory,
    "compartment_coordination_method_diagram.png"
  ),
  plot = coordination_method_figure,
  width = 13,
  height = 10,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_directory,
    "compartment_coordination_method_diagram.pdf"
  ),
  plot = coordination_method_figure,
  width = 13,
  height = 10,
  units = "in",
  bg = "white"
)
