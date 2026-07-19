# 1) Prepare microbiome data and pairwise distances ####

library(tidyverse)
library(vegan)

# 0) File paths ####

data_directory <- "../Data/"
output_directory <- file.path("../Output/", "Poster_Common_5_Sites")
intermediate_directory <- file.path(output_directory, "Intermediate")

dir.create(output_directory, showWarnings = FALSE, recursive = TRUE)
dir.create(intermediate_directory, showWarnings = FALSE, recursive = TRUE)

metadata_file <- file.path(
  data_directory,
  "FPCP_metadata.csv"
)

taxa_16S_file <- file.path(
  data_directory,
  "16S_emu-combined-tax_id.tsv"
)

taxa_799F_file <- file.path(
  data_directory,
  "799F_emu-combined-tax_id.tsv"
)

# 1) Prepare genus-level abundance and metadata ####

## 1a) Read metadata and EMU taxonomy tables ####

metadata <- read.csv(
  metadata_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

taxa16S <- read.delim(
  taxa_16S_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

taxa799F <- read.delim(
  taxa_799F_file,
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

## 1b) Standardize taxonomic identifiers ####

taxa16S$tax_id <- as.character(taxa16S$tax_id)
taxa799F$tax_id <- as.character(taxa799F$tax_id)

## 1c) Join the 16S and 799F abundance tables ####

taxa <- full_join(
  taxa16S,
  taxa799F,
  by = "tax_id",
  suffix = c(".16S", ".799F")
)

## 1d) Reconcile duplicated taxonomy annotations ####

taxonomy_ranks <- c(
  "species",
  "genus",
  "family",
  "order",
  "class",
  "phylum",
  "superkingdom"
)

for (rank in taxonomy_ranks) {
  rank_16S <- paste0(rank, ".16S")
  rank_799F <- paste0(rank, ".799F")

  taxa[[rank]] <- coalesce(
    taxa[[rank_16S]],
    taxa[[rank_799F]]
  )
}

## 1e) Collapse abundances to genus ####

taxa_g <- taxa %>%
  mutate(
    genus = case_when(
      is.na(genus) | genus == "" ~
        paste0("Unclassified_", tax_id),
      TRUE ~ genus
    )
  ) %>%
  group_by(genus) %>%
  summarise(
    across(
      matches("\\.filtered$"),
      ~ sum(as.numeric(.x), na.rm = TRUE)
    ),
    .groups = "drop"
  )

taxa_f <- taxa %>%
  mutate(
    family = case_when(
      is.na(family) | family == "" ~
        paste0("Unclassified_", tax_id),
      TRUE ~ family
    )
  ) %>%
  group_by(family) %>%
  summarise(
    across(
      matches("\\.filtered$"),
      ~ sum(as.numeric(.x), na.rm = TRUE)
    ),
    .groups = "drop"
  )

## 1f) Create the sample-by-genus matrix ####

sample_columns <- grep(
  "\\.filtered$",
  colnames(taxa_g),
  value = TRUE
)

if (length(sample_columns) == 0) {
  stop(
    "No sample columns ending in '.filtered' were found."
  )
}

taxa_abundance <- taxa_g %>%
  select(all_of(sample_columns))

taxa_abundance[is.na(taxa_abundance)] <- 0
taxa_abundance[] <- lapply(
  taxa_abundance,
  as.numeric
)

taxa_abundance_matrix <- t(
  as.matrix(taxa_abundance)
)

rownames(taxa_abundance_matrix) <- sub(
  "\\.filtered$",
  "",
  rownames(taxa_abundance_matrix)
)

colnames(taxa_abundance_matrix) <- taxa_g$genus
storage.mode(taxa_abundance_matrix) <- "numeric"

## 1g) Match metadata and abundance samples ####

metadata$Sample_Index <- as.character(
  metadata$Sample_Index
)

common_samples <- intersect(
  metadata$Sample_Index,
  rownames(taxa_abundance_matrix)
)

if (length(common_samples) == 0) {
  stop(
    paste(
      "No sample IDs were shared by metadata",
      "and the abundance matrix."
    )
  )
}

metadata_matched <- metadata %>%
  filter(
    Sample_Index %in% common_samples
  ) %>%
  arrange(
    match(Sample_Index, common_samples)
  )

rownames(metadata_matched) <-
  metadata_matched$Sample_Index

taxa_abundance_matched <-
  taxa_abundance_matrix[
    rownames(metadata_matched),
    ,
    drop = FALSE
  ]

stopifnot(
  identical(
    rownames(metadata_matched),
    rownames(taxa_abundance_matched)
  )
)

## 1h) Assign microbiome compartments ####

metadata_matched <- metadata_matched %>%
  mutate(
    Site = as.character(Sample_Site),
    Plant = as.character(Sample_Plant),

    # Plant numbers repeat among sites, so site must be included.
    Plant_ID = paste(
      Site,
      Plant,
      sep = "_"
    ),

    Compartment = case_when(
      Sample_Type1 == "Leaf" &
        Sample_Type2 == "Endophyte" ~
        "Leaf_EN",

      Sample_Type1 == "Leaf" &
        Sample_Type2 == "Epiphyte" ~
        "Leaf_EP",

      Sample_Type1 == "Fruit" &
        Sample_Type2 == "Endophyte" ~
        "Fruit_EN",

      Sample_Type1 == "Fruit" &
        Sample_Type2 == "Epiphyte" ~
        "Fruit_EP",

      Sample_Type1 == "Root" ~
        "Root",

      Sample_Type1 == "Soil" &
        Sample_Type2 == "Plant" ~
        "Soil_p",

      Sample_Type1 == "Soil" &
        Sample_Type2 == "Meter" ~
        "Soil_1m",

      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(Compartment)
  )

## This analysis assumes one sample for each plant-compartment combination.
duplicate_compartments <- metadata_matched %>%
  count(
    Site,
    Plant,
    Plant_ID,
    Compartment,
    name = "N"
  ) %>%
  filter(
    N > 1
  )

if (nrow(duplicate_compartments) > 0) {
  stop(
    paste(
      "At least one plant-compartment combination contains",
      "multiple samples. Resolve or aggregate those replicates",
      "before calculating plant-level Bray-Curtis distances."
    )
  )
}

## 1i) Define compartment order and labels ####

compartment_order <- c(
  "Soil_1m",
  "Soil_p",
  "Root",
  "Leaf_EN",
  "Leaf_EP",
  "Fruit_EN",
  "Fruit_EP"
)


## 1j) Restrict the primary analysis to complete sites ####
##
## This poster-stage analysis uses only sites represented in every
## compartment so that soil, rhizosphere, root, leaf, and fruit results
## are based on the same set of sites.
##
## When fruit processing is complete for all eight sites, remove or
## disable this block to restore the full analysis.

site_compartment_coverage <- metadata_matched %>%
  distinct(
    Site,
    Compartment
  ) %>%
  count(
    Site,
    name = "N_Compartments"
  ) %>%
  arrange(
    Site
  )

complete_sites <- site_compartment_coverage %>%
  filter(
    N_Compartments == length(compartment_order)
  ) %>%
  pull(
    Site
  ) %>%
  as.character()

if (length(complete_sites) < 2) {
  stop(
    paste(
      "Fewer than two sites contain all seven compartments.",
      "Inspect site_compartment_coverage."
    )
  )
}

metadata_matched <- metadata_matched %>%
  filter(
    Site %in% complete_sites
  ) %>%
  arrange(
    match(
      Sample_Index,
      rownames(taxa_abundance_matched)
    )
  )

taxa_abundance_matched <- taxa_abundance_matched[
  metadata_matched$Sample_Index,
  ,
  drop = FALSE
]

stopifnot(
  identical(
    metadata_matched$Sample_Index,
    rownames(taxa_abundance_matched)
  )
)

# Confirm retained sampling after restriction. ###
retained_sample_summary <- metadata_matched %>%
  count(
    Site,
    Compartment,
    name = "N_Samples"
  ) %>%
  arrange(
    Site,
    factor(
      Compartment,
      levels = compartment_order
    )
  )

write_csv(
  site_compartment_coverage,
  file.path(
    output_directory,
    "site_compartment_coverage_before_restriction.csv"
  )
)

write_csv(
  retained_sample_summary,
  file.path(
    output_directory,
    "retained_samples_complete_sites.csv"
  )
)

## 2a) Define all-pair Bray-Curtis calculation ####

calculate_all_pair_bray <- function(
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
      as.numeric(Plant),
      Plant_ID
    )

  if (nrow(compartment_metadata) < 2) {
    return(
      tibble()
    )
  }

  compartment_abundance <- abundance_matrix[
    compartment_metadata$Sample_Index,
    ,
    drop = FALSE
  ]

  # Remove genera absent from every sample in this compartment.
  compartment_abundance <-
    compartment_abundance[
      ,
      colSums(compartment_abundance) > 0,
      drop = FALSE
    ]

  if (ncol(compartment_abundance) == 0) {
    return(
      tibble()
    )
  }

  compartment_bray <- as.matrix(
    vegdist(
      compartment_abundance,
      method = "bray"
    )
  )

  plant_pairs <- combn(
    seq_len(nrow(compartment_metadata)),
    2
  )

  tibble(
    Site_1 =
      compartment_metadata$Site[
        plant_pairs[1, ]
      ],

    Site_2 =
      compartment_metadata$Site[
        plant_pairs[2, ]
      ],

    Plant_1 =
      compartment_metadata$Plant[
        plant_pairs[1, ]
      ],

    Plant_2 =
      compartment_metadata$Plant[
        plant_pairs[2, ]
      ],

    Plant_1_ID =
      compartment_metadata$Plant_ID[
        plant_pairs[1, ]
      ],

    Plant_2_ID =
      compartment_metadata$Plant_ID[
        plant_pairs[2, ]
      ],

    Bray_Curtis =
      compartment_bray[
        cbind(
          plant_pairs[1, ],
          plant_pairs[2, ]
        )
      ],

    Compartment = compartment_name
  ) %>%
    mutate(
      Same_Site = Site_1 == Site_2
    )
}

## 2b) Calculate compartment-specific distances ####

bray_long <- map_dfr(
  compartment_order,
  calculate_all_pair_bray,
  metadata_table = metadata_matched,
  abundance_matrix = taxa_abundance_matched
)

## 2c) Construct the complete plant-pair index ####

plant_lookup <- metadata_matched %>%
  distinct(
    Site,
    Plant,
    Plant_ID
  ) %>%
  arrange(
    Site,
    as.numeric(Plant),
    Plant_ID
  )

if (nrow(plant_lookup) < 2) {
  stop(
    "Fewer than two unique plants were available."
  )
}

all_pair_index <- combn(
  seq_len(nrow(plant_lookup)),
  2
)

all_pairs <- tibble(
  Site_1 =
    plant_lookup$Site[
      all_pair_index[1, ]
    ],

  Site_2 =
    plant_lookup$Site[
      all_pair_index[2, ]
    ],

  Plant_1 =
    plant_lookup$Plant[
      all_pair_index[1, ]
    ],

  Plant_2 =
    plant_lookup$Plant[
      all_pair_index[2, ]
    ],

  Plant_1_ID =
    plant_lookup$Plant_ID[
      all_pair_index[1, ]
    ],

  Plant_2_ID =
    plant_lookup$Plant_ID[
      all_pair_index[2, ]
    ]
) %>%
  mutate(
    Same_Site = Site_1 == Site_2,
    Pair = paste(
      Plant_1_ID,
      "vs",
      Plant_2_ID
    )
  )

## 2d) Assemble the plant-pair-by-compartment table ####

bray_pair_table <- all_pairs %>%
  left_join(
    bray_long,
    by = c(
      "Site_1",
      "Site_2",
      "Plant_1",
      "Plant_2",
      "Plant_1_ID",
      "Plant_2_ID",
      "Same_Site"
    )
  ) %>%
  select(
    Site_1,
    Site_2,
    Same_Site,
    Plant_1,
    Plant_2,
    Plant_1_ID,
    Plant_2_ID,
    Pair,
    Compartment,
    Bray_Curtis
  ) %>%
  pivot_wider(
    names_from = Compartment,
    values_from = Bray_Curtis
  ) %>%
  select(
    Site_1,
    Site_2,
    Same_Site,
    Plant_1,
    Plant_2,
    Plant_1_ID,
    Plant_2_ID,
    Pair,
    all_of(compartment_order)
  )

write_csv(
  bray_pair_table,
  file.path(
    output_directory,
    "all_plant_pair_bray_table.csv"
  )
)

# 3) Save cleaned analysis objects ####

clean_data_objects <- list(
  metadata_matched = metadata_matched,
  taxa_abundance_matched = taxa_abundance_matched,
  bray_long = bray_long,
  bray_pair_table = bray_pair_table,
  compartment_order = compartment_order
)

saveRDS(
  clean_data_objects,
  file.path(intermediate_directory, "clean_data_objects.rds")
)

