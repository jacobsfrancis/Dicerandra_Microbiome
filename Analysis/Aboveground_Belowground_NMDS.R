# NMDS ordinations for aboveground and belowground microbiomes ####
#
# This script uses the cleaned objects created by 01_Data_Cleaning.R.
#
# Aboveground compartments:
#   Leaf endophytes, leaf epiphytes, fruit endophytes, fruit epiphytes
#
# Belowground compartments:
#   Bulk soil, rhizosphere, roots
#
# The two groups are ordinated separately because they were generated using
# different sequencing approaches.

library(tidyverse)
library(vegan)
library(patchwork)

# 0) File paths ####

output_directory <- file.path(
  "../Output/",
  "Poster_Common_5_Sites"
)

intermediate_directory <- file.path(
  output_directory,
  "Intermediate"
)

clean_data_file <- file.path(
  intermediate_directory,
  "clean_data_objects.rds"
)

if (!file.exists(clean_data_file)) {
  stop(
    paste(
      "Could not find:",
      clean_data_file,
      "\nRun 01_Data_Cleaning.R first."
    )
  )
}

# 1) Load cleaned data ####

clean_data_objects <- readRDS(
  clean_data_file
)

metadata_matched <- clean_data_objects$metadata_matched
taxa_abundance_matched <- clean_data_objects$taxa_abundance_matched

stopifnot(
  identical(
    metadata_matched$Sample_Index,
    rownames(taxa_abundance_matched)
  )
)

# 2) Define compartment order, labels, colors, and shapes ####

aboveground_compartments <- c(
  "Leaf_EN",
  "Leaf_EP",
  "Fruit_EN",
  "Fruit_EP"
)

belowground_compartments <- c(
  "Soil_1m",
  "Soil_p",
  "Root"
)

compartment_labels <- c(
  Soil_1m  = "Bulk soil",
  Soil_p   = "Rhizosphere",
  Root     = "Root",
  Leaf_EN  = "Leaf endophyte",
  Leaf_EP  = "Leaf epiphyte",
  Fruit_EN = "Fruit endophyte",
  Fruit_EP = "Fruit epiphyte"
)

compartment_colors <- c(
  Soil_1m  = "#C83A06", # bulk soil: red
  Soil_p   = "#DD6038", # rhizosphere: orange
  Root     = "#F2C54D", # root: yellow
  Leaf_EN  = "#476C35", # leaf endophyte: dark green
  Leaf_EP  = "#B6CD78", # leaf epiphyte: light green
  Fruit_EN = "#355EAA", # fruit endophyte: dark blue
  Fruit_EP = "#6EA4CE"  # fruit epiphyte: light blue
)

compartment_shapes <- c(
  Soil_1m  = 22,
  Soil_p   = 22,
  Root     = 23,
  Leaf_EN  = 24,
  Leaf_EP  = 24,
  Fruit_EN = 21,
  Fruit_EP = 21
)

# 3) Set the poster plotting theme ####

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
      legend.text = element_text(
        size = 16
      ),
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

# 4) Function to run one NMDS ####

run_compartment_nmds <- function(
    compartment_subset,
    metadata_table,
    abundance_matrix,
    random_seed = 20260718
) {

  nmds_metadata <- metadata_table %>%
    filter(
      Compartment %in% compartment_subset
    ) %>%
    mutate(
      Compartment = factor(
        Compartment,
        levels = compartment_subset
      )
    )

  nmds_abundance <- abundance_matrix[
    nmds_metadata$Sample_Index,
    ,
    drop = FALSE
  ]

  keep_samples <- rowSums(
    nmds_abundance,
    na.rm = TRUE
  ) > 0

  nmds_metadata <- nmds_metadata[
    keep_samples,
    ,
    drop = FALSE
  ]

  nmds_abundance <- nmds_abundance[
    keep_samples,
    ,
    drop = FALSE
  ]

  keep_taxa <- colSums(
    nmds_abundance,
    na.rm = TRUE
  ) > 0

  nmds_abundance <- nmds_abundance[
    ,
    keep_taxa,
    drop = FALSE
  ]

  if (nrow(nmds_abundance) < 4) {
    stop(
      "Fewer than four nonempty samples remained for one NMDS."
    )
  }

  if (ncol(nmds_abundance) < 2) {
    stop(
      "Fewer than two taxa remained for one NMDS."
    )
  }

  set.seed(
    random_seed
  )

  nmds_fit <- vegan::metaMDS(
    nmds_abundance,
    distance = "bray",
    k = 2,
    trymax = 200,
    autotransform = FALSE,
    trace = FALSE
  )

  nmds_scores <- vegan::scores(
    nmds_fit,
    display = "sites"
  ) %>%
    as.data.frame() %>%
    rownames_to_column(
      "Sample_Index"
    ) %>%
    as_tibble() %>%
    left_join(
      nmds_metadata %>%
        transmute(
          Sample_Index = as.character(Sample_Index),
          Site,
          Plant,
          Plant_ID,
          Compartment
        ),
      by = "Sample_Index"
    )

  list(
    fit = nmds_fit,
    scores = nmds_scores
  )
}

# 5) Run separate aboveground and belowground NMDS ordinations ####

aboveground_nmds <- run_compartment_nmds(
  compartment_subset = aboveground_compartments,
  metadata_table = metadata_matched,
  abundance_matrix = taxa_abundance_matched,
  random_seed = 20260718
)

belowground_nmds <- run_compartment_nmds(
  compartment_subset = belowground_compartments,
  metadata_table = metadata_matched,
  abundance_matrix = taxa_abundance_matched,
  random_seed = 20260719
)

# 6) Function to draw one NMDS panel ####

make_nmds_plot <- function(
    nmds_object,
    panel_title,
    compartment_subset
) {

  score_data <- nmds_object$scores %>%
    mutate(
      Compartment = factor(
        Compartment,
        levels = compartment_subset
      )
    )

  stress_label <- sprintf(
    "Stress = %.3f",
    nmds_object$fit$stress
  )

  ggplot(
    score_data,
    aes(
      x = NMDS1,
      y = NMDS2,
      fill = Compartment,
      color = Compartment,
      shape = Compartment
    )
  ) +

    stat_ellipse(
      aes(
        group = Compartment
      ),
      type = "t",
      level = 0.95,
      linewidth = 1.1,
      alpha = 0.85,
      show.legend = FALSE
    ) +

    geom_point(
      size = 3.4,
      stroke = 0.8,
      alpha = 0.72
    ) +

    annotate(
      geom = "text",
      x = Inf,
      y = Inf,
      label = stress_label,
      hjust = 1.05,
      vjust = 1.25,
      size = 16 / .pt,
      family = "Helvetica",
      color = "#262626"
    ) +

    scale_fill_manual(
      values = compartment_colors[compartment_subset],
      labels = compartment_labels[compartment_subset],
      drop = FALSE
    ) +

    scale_color_manual(
      values = compartment_colors[compartment_subset],
      labels = compartment_labels[compartment_subset],
      drop = FALSE
    ) +

    scale_shape_manual(
      values = compartment_shapes[compartment_subset],
      labels = compartment_labels[compartment_subset],
      drop = FALSE
    ) +

    labs(
      title = panel_title,
      x = "NMDS1",
      y = "NMDS2"
    ) +

    guides(
      fill = guide_legend(
        override.aes = list(
          alpha = 1,
          size = 4
        )
      ),
      color = "none",
      shape = "none"
    ) +

    theme(
      legend.position = "right"
    )
}

# 7) Draw the two ordinations ####

aboveground_plot <- make_nmds_plot(
  nmds_object = aboveground_nmds,
  panel_title = "A. Aboveground",
  compartment_subset = aboveground_compartments
)

belowground_plot <- make_nmds_plot(
  nmds_object = belowground_nmds,
  panel_title = "B. Belowground",
  compartment_subset = belowground_compartments
)

# 8) Stack and save the figure ####

community_nmds_figure <-
  aboveground_plot +
  belowground_plot +
  patchwork::plot_layout(
    nrow=1,
    guides = "collect",
    heights = c(
      1,
      1
    )
  ) &
  theme(
    legend.position = "none"
  )

community_nmds_figure

ggsave(
  filename = file.path(
    output_directory,
    "aboveground_belowground_NMDS.png"),
  plot = community_nmds_figure,
  width = 9.75,
  height = 4,
  units = "in",
  dpi = 600,
  bg = "transparent"
)

ggsave(
  filename = file.path(
    output_directory,
    "aboveground_belowground_NMDS.pdf"
  ),
  plot = community_nmds_figure,
  width = 9.75,
  height = 4,
  units = "in",
  bg = "transparent"
)

# 9) Save NMDS scores and stress values ####

write_csv(
  aboveground_nmds$scores,
  file.path(
    output_directory,
    "aboveground_NMDS_scores.csv"
  )
)

write_csv(
  belowground_nmds$scores,
  file.path(
    output_directory,
    "belowground_NMDS_scores.csv"
  )
)

nmds_stress_summary <- tibble(
  Ordination = c(
    "Aboveground",
    "Belowground"
  ),
  Stress = c(
    aboveground_nmds$fit$stress,
    belowground_nmds$fit$stress
  )
)

write_csv(
  nmds_stress_summary,
  file.path(
    output_directory,
    "aboveground_belowground_NMDS_stress.csv"
  )
)

#Permanova analysis####

## PERMANOVA: Aboveground compartments ####

above_metadata <- metadata_matched %>%
  filter(
    Compartment %in% aboveground_compartments
  )

above_abundance <- taxa_abundance_matched[
  above_metadata$Sample_Index,
  ,
  drop = FALSE
]

above_bray <- vegdist(
  above_abundance,
  method = "bray"
)

above_permanova <- adonis2(
  above_bray ~ Compartment,
  data = above_metadata,
  permutations = 9999
)

above_permanova


## PERMANOVA: Belowground compartments ####

below_metadata <- metadata_matched %>%
  filter(
    Compartment %in% belowground_compartments
  )

below_abundance <- taxa_abundance_matched[
  below_metadata$Sample_Index,
  ,
  drop = FALSE
]

below_bray <- vegdist(
  below_abundance,
  method = "bray"
)

below_permanova <- adonis2(
  below_bray ~ Compartment,
  data = below_metadata,
  permutations = 9999
)

below_permanova


#Pairwise
pairwise_permanova <- function(metadata,
                               abundance,
                               group_column,
                               permutations = 9999){
  
  groups <- sort(unique(metadata[[group_column]]))
  
  comparisons <- combn(groups, 2, simplify = FALSE)
  
  results <- purrr::map_dfr(comparisons, function(comp){
    
    keep <- metadata[[group_column]] %in% comp
    
    md <- metadata[keep, ]
    
    abund <- abundance[
      md$Sample_Index,
      ,
      drop = FALSE
    ]
    
    bray <- vegdist(
      abund,
      method = "bray"
    )
    
    fit <- adonis2(
      bray ~ md[[group_column]],
      permutations = permutations
    )
    
    tibble(
      Group1 = comp[1],
      Group2 = comp[2],
      F = fit$F[1],
      R2 = fit$R2[1],
      P = fit$`Pr(>F)`[1]
    )
    
  })
  
  results$P_adj <-
    p.adjust(
      results$P,
      method = "BH"
    )
  
  results
}

above_pairwise <-
  pairwise_permanova(
    metadata = above_metadata,
    abundance = above_abundance,
    group_column = "Compartment"
  )

below_pairwise <-
  pairwise_permanova(
    metadata = below_metadata,
    abundance = below_abundance,
    group_column = "Compartment"
  )
