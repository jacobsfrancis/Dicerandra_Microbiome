# 8) Plot observed genus richness ####

library(tidyverse)

# 0) File paths ####

data_directory <- "../Data/"
output_directory <- file.path("../Output/", "Poster_Common_5_Sites")
intermediate_directory <- file.path(output_directory, "Intermediate")

dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
dir.create(intermediate_directory, showWarnings = FALSE, recursive = TRUE)

# 1) Load cleaned data ####

clean_data_objects <- readRDS(
  file.path(intermediate_directory, "clean_data_objects.rds")
)
list2env(clean_data_objects, envir = .GlobalEnv)

# 2) Calculate sample read depth and observed richness ####

sample_richness <- tibble(
  Sample_Index = rownames(
    taxa_abundance_matched
  ),

  Richness = rowSums(
    taxa_abundance_matched > 0
  ),

  Total_Reads = rowSums(
    taxa_abundance_matched
  )
) %>%
  left_join(
    metadata_matched %>%
      select(
        Sample_Index,
        Site,
        Plant_ID,
        Compartment
      ),
    by = "Sample_Index"
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
    )
  )

write_csv(
  sample_richness,
  file.path(
    output_directory,
    "sample_read_depth_and_richness.csv"
  )
)


# 3) Draw and save observed-richness plot ####

richness_plot <- ggplot(
  sample_richness,
  aes(
    x = Richness,
    y = Compartment_Label
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.alpha = 0.7
  ) +
  labs(
    title = "Observed genus richness across microbiome compartments",
    x = "Observed genera",
    y = NULL
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    plot.title = element_text(
      face = "bold"
    ),
    axis.text.y = element_text(
      face = "bold"
    )
  )

ggsave(
  filename = file.path(
    output_directory,
    "observed_genus_richness_by_compartment.png"
  ),
  plot = richness_plot,
  width = 8,
  height = 5.5,
  units = "in",
  dpi = 400,
  bg = "white"
)

