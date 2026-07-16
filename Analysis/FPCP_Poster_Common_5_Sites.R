## ============================================================
## COORDINATION OF MICROBIOME DISSIMILARITY ACROSS COMPARTMENTS####
## ALL-PAIR SENSITIVITY ANALYSIS####
## ============================================================
##
## Biological question:
## Across all sampled plants, are pairs of plants that differ more strongly
## in one microbiome compartment also more different in another compartment,
## after accounting for whether the two plants came from the same site?
##
## Analysis:
## 1. Combine the 16S and 799F EMU abundance tables.
## 2. Collapse abundance to genus.
## 3. Calculate Bray-Curtis dissimilarities among ALL plant pairs,
##    separately for each compartment.
## 4. Fit one standardized mixed model for each unique compartment pair:
##
##      Response dissimilarity ~ Predictor dissimilarity + Same_Site +
##        (1 | Site_1) + (1 | Site_2) +
##        (1 | Plant_1_ID) + (1 | Plant_2_ID)
##
##    Same_Site controls the average difference between within-site and
##    between-site pairs. The crossed site and plant random intercepts account
##    for repeated use of the same sites and plants across pairwise comparisons.
##
## 5. Adjust the 21 predictor p-values using Benjamini-Hochberg FDR correction.
## 6. Produce:
##    - a color-coded matrix of standardized predictor coefficients
##    - a coefficient plot with estimates and 95% confidence intervals
##    - a separate table containing the Same_Site fixed effect from each model
##
## IMPORTANT INTERPRETATION:
## The plotted coefficient for Predictor_z estimates the association between
## dissimilarities across the two compartments after controlling for whether
## the pair is within the same site and accounting for repeated use of sites
## and plants.
##
## The heatmap contains standardized mixed-model coefficients, not raw
## Pearson correlations.
## ============================================================

library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(broom.mixed)

## ------------------------------------------------------------
## 0. File paths and output directory####
## ------------------------------------------------------------

data_directory <- "../Data/"
output_directory <- file.path(
  "../Output/",
  "Poster_All_Available_Plants"
)

dir.create(
  output_directory,
  showWarnings = FALSE,
  recursive = TRUE
)

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

## ------------------------------------------------------------
## 1. Read metadata and both EMU taxonomy tables
## ------------------------------------------------------------

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

## ------------------------------------------------------------
## 2. Make tax_id the same type in both files
## ------------------------------------------------------------

taxa16S$tax_id <- as.character(taxa16S$tax_id)
taxa799F$tax_id <- as.character(taxa799F$tax_id)

## ------------------------------------------------------------
## 3. Combine full-length 16S and 799F tables by tax_id
## ------------------------------------------------------------

taxa <- full_join(
  taxa16S,
  taxa799F,
  by = "tax_id",
  suffix = c(".16S", ".799F")
)

## ------------------------------------------------------------
## 4. Coalesce duplicated taxonomy-rank columns
## ------------------------------------------------------------

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

## ------------------------------------------------------------
## 5. Collapse abundance values to genus
## ------------------------------------------------------------

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

## ------------------------------------------------------------
## 6. Create samples x genera abundance matrix
## ------------------------------------------------------------

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

## ------------------------------------------------------------
## 7. Match metadata to the abundance matrix
## ------------------------------------------------------------

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

## ------------------------------------------------------------
## 8. Assign samples to analysis compartments
## ------------------------------------------------------------

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
  print(duplicate_compartments)

  stop(
    paste(
      "At least one plant-compartment combination contains",
      "multiple samples. Resolve or aggregate those replicates",
      "before calculating plant-level Bray-Curtis distances."
    )
  )
}

## ------------------------------------------------------------
## 9. Compartment order
## ------------------------------------------------------------

compartment_order <- c(
  "Soil_1m",
  "Soil_p",
  "Root",
  "Leaf_EN",
  "Leaf_EP",
  "Fruit_EN",
  "Fruit_EP"
)


## ------------------------------------------------------------
## 9A. Restrict all analyses to sites with all seven compartments
## ------------------------------------------------------------
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

print(
  site_compartment_coverage
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

message(
  "Restricting analyses to complete sites: ",
  paste(
    complete_sites,
    collapse = ", "
  )
)

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

print(
  retained_sample_summary
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

## ------------------------------------------------------------
## 10. Function: all-pair Bray-Curtis for one compartment
## ------------------------------------------------------------

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

## ------------------------------------------------------------
## 11. Calculate all-pair distances for every compartment
## ------------------------------------------------------------

bray_long <- map_dfr(
  compartment_order,
  calculate_all_pair_bray,
  metadata_table = metadata_matched,
  abundance_matrix = taxa_abundance_matched
)

## ------------------------------------------------------------
## 12. Construct every expected pair among sampled plants
## ------------------------------------------------------------

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

## ------------------------------------------------------------
## 13. Create final all-pair x compartment table
## ------------------------------------------------------------

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

print(
  bray_pair_table %>%
    count(
      Same_Site
    )
)

## ============================================================
## STANDARDIZED DYADIC LINEAR MIXED MODELS
## ============================================================

## ------------------------------------------------------------
## 14. Example model: Leaf endophyte versus root dissimilarity
## ------------------------------------------------------------

leaf_en_root_data <- bray_pair_table %>%
  select(
    Site_1,
    Site_2,
    Same_Site,
    Plant_1_ID,
    Plant_2_ID,
    Leaf_EN,
    Root
  ) %>%
  filter(
    is.finite(Leaf_EN),
    is.finite(Root)
  ) %>%
  mutate(
    Site_1 = factor(Site_1),
    Site_2 = factor(Site_2),

    Same_Site = factor(
      Same_Site,
      levels = c(FALSE, TRUE),
      labels = c(
        "Different sites",
        "Same site"
      )
    ),

    Plant_1_ID = factor(Plant_1_ID),
    Plant_2_ID = factor(Plant_2_ID),

    Leaf_EN_z = as.numeric(
      scale(Leaf_EN)
    ),

    Root_z = as.numeric(
      scale(Root)
    )
  )

leaf_en_root_model <- lmer(
  Leaf_EN_z ~
    Root_z +
    Same_Site +
    (1 | Site_1) +
    (1 | Site_2) +
    (1 | Plant_1_ID) +
    (1 | Plant_2_ID),

  data = leaf_en_root_data,
  REML = FALSE,

  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(
      maxfun = 200000
    )
  )
)

print(
  summary(leaf_en_root_model)
)

## ------------------------------------------------------------
## 15. Function to fit one compartment-pair model
## ------------------------------------------------------------

fit_distance_lmm <- function(
    response_name,
    predictor_name,
    data
) {
  model_data <- data %>%
    select(
      Site_1,
      Site_2,
      Same_Site,
      Plant_1_ID,
      Plant_2_ID,
      all_of(response_name),
      all_of(predictor_name)
    ) %>%
    rename(
      Response_distance =
        all_of(response_name),

      Predictor_distance =
        all_of(predictor_name)
    ) %>%
    filter(
      is.finite(Response_distance),
      is.finite(Predictor_distance)
    ) %>%
    mutate(
      Site_1 = factor(Site_1),
      Site_2 = factor(Site_2),

      Same_Site = factor(
        Same_Site,
        levels = c(FALSE, TRUE),
        labels = c(
          "Different sites",
          "Same site"
        )
      ),

      Plant_1_ID = factor(Plant_1_ID),
      Plant_2_ID = factor(Plant_2_ID)
    )

  number_of_pairs <- nrow(model_data)

  number_of_within_site_pairs <- sum(
    model_data$Same_Site == "Same site"
  )

  number_of_between_site_pairs <- sum(
    model_data$Same_Site == "Different sites"
  )

  number_of_sites <- n_distinct(
    c(
      as.character(model_data$Site_1),
      as.character(model_data$Site_2)
    )
  )

  number_of_plants <- n_distinct(
    c(
      as.character(model_data$Plant_1_ID),
      as.character(model_data$Plant_2_ID)
    )
  )

  response_sd <- sd(
    model_data$Response_distance
  )

  predictor_sd <- sd(
    model_data$Predictor_distance
  )

  empty_result <- function(
      status_message
  ) {
    tibble(
      Response = response_name,
      Predictor = predictor_name,

      Estimate = NA_real_,
      Standard_Error = NA_real_,
      Degrees_of_Freedom = NA_real_,
      Test_Statistic = NA_real_,
      P_Value = NA_real_,
      CI_Lower = NA_real_,
      CI_Upper = NA_real_,

      Same_Site_Estimate = NA_real_,
      Same_Site_Standard_Error = NA_real_,
      Same_Site_Degrees_of_Freedom = NA_real_,
      Same_Site_Test_Statistic = NA_real_,
      Same_Site_P_Value = NA_real_,
      Same_Site_CI_Lower = NA_real_,
      Same_Site_CI_Upper = NA_real_,

      N_Pairs = number_of_pairs,

      N_Within_Site_Pairs =
        number_of_within_site_pairs,

      N_Between_Site_Pairs =
        number_of_between_site_pairs,

      N_Sites = number_of_sites,
      N_Plants = number_of_plants,

      Site_1_Variance = NA_real_,
      Site_2_Variance = NA_real_,
      Plant_1_Variance = NA_real_,
      Plant_2_Variance = NA_real_,
      Residual_Variance = NA_real_,

      Singular = NA,
      Convergence_Message = NA_character_,
      Model_Status = status_message
    )
  }

  if (
    number_of_pairs < 5 ||
    number_of_sites < 2 ||
    number_of_plants < 3 ||
    n_distinct(model_data$Same_Site) < 2 ||
    is.na(response_sd) ||
    is.na(predictor_sd) ||
    response_sd == 0 ||
    predictor_sd == 0
  ) {
    return(
      empty_result(
        "Insufficient data or variation"
      )
    )
  }

  model_data <- model_data %>%
    mutate(
      Response_z = as.numeric(
        scale(Response_distance)
      ),

      Predictor_z = as.numeric(
        scale(Predictor_distance)
      )
    )

  model_fit <- tryCatch(
    lmer(
      Response_z ~
        Predictor_z +
        Same_Site +
        (1 | Site_1) +
        (1 | Site_2) +
        (1 | Plant_1_ID) +
        (1 | Plant_2_ID),

      data = model_data,
      REML = FALSE,

      control = lmerControl(
        optimizer = "bobyqa",
        optCtrl = list(
          maxfun = 200000
        )
      )
    ),

    error = function(e) e
  )

  if (inherits(model_fit, "error")) {
    return(
      empty_result(
        model_fit$message
      )
    )
  }

  fixed_effects <- broom.mixed::tidy(
    model_fit,
    effects = "fixed",
    conf.int = TRUE,
    conf.level = 0.95
  )

  predictor_effect <- fixed_effects %>%
    filter(
      term == "Predictor_z"
    )

  same_site_effect <- fixed_effects %>%
    filter(
      str_detect(
        term,
        "^Same_Site"
      )
    )

  if (nrow(predictor_effect) != 1) {
    return(
      empty_result(
        paste(
          "Predictor coefficient",
          "was not returned"
        )
      )
    )
  }

  if (nrow(same_site_effect) != 1) {
    return(
      empty_result(
        paste(
          "Same_Site coefficient",
          "was not returned"
        )
      )
    )
  }

  variance_table <- as.data.frame(
    VarCorr(model_fit)
  )

  get_variance <- function(
      group_name
  ) {
    value <- variance_table %>%
      filter(
        grp == group_name
      ) %>%
      pull(
        vcov
      )

    if (length(value) == 0) {
      NA_real_
    } else {
      value[1]
    }
  }

  convergence_messages <-
    model_fit@optinfo$conv$lme4$messages

  convergence_message <-
    if (is.null(convergence_messages)) {
      NA_character_
    } else {
      paste(
        convergence_messages,
        collapse = "; "
      )
    }

  tibble(
    Response = response_name,
    Predictor = predictor_name,

    Estimate =
      predictor_effect$estimate,

    Standard_Error =
      predictor_effect$std.error,

    Degrees_of_Freedom =
      predictor_effect$df,

    Test_Statistic =
      predictor_effect$statistic,

    P_Value =
      predictor_effect$p.value,

    CI_Lower =
      predictor_effect$conf.low,

    CI_Upper =
      predictor_effect$conf.high,

    Same_Site_Estimate =
      same_site_effect$estimate,

    Same_Site_Standard_Error =
      same_site_effect$std.error,

    Same_Site_Degrees_of_Freedom =
      same_site_effect$df,

    Same_Site_Test_Statistic =
      same_site_effect$statistic,

    Same_Site_P_Value =
      same_site_effect$p.value,

    Same_Site_CI_Lower =
      same_site_effect$conf.low,

    Same_Site_CI_Upper =
      same_site_effect$conf.high,

    N_Pairs = number_of_pairs,

    N_Within_Site_Pairs =
      number_of_within_site_pairs,

    N_Between_Site_Pairs =
      number_of_between_site_pairs,

    N_Sites = number_of_sites,
    N_Plants = number_of_plants,

    Site_1_Variance =
      get_variance("Site_1"),

    Site_2_Variance =
      get_variance("Site_2"),

    Plant_1_Variance =
      get_variance("Plant_1_ID"),

    Plant_2_Variance =
      get_variance("Plant_2_ID"),

    Residual_Variance =
      get_variance("Residual"),

    Singular =
      isSingular(
        model_fit,
        tol = 1e-4
      ),

    Convergence_Message =
      convergence_message,

    Model_Status = "Success"
  )
}

## ------------------------------------------------------------
## 16. Create all 21 unique compartment comparisons
## ------------------------------------------------------------
##
## The earlier compartment is used as the predictor and the later
## compartment as the response. Each unordered pair appears once.

model_grid <- combn(
  compartment_order,
  2,
  simplify = FALSE
) %>%
  map_dfr(
    ~ tibble(
      Predictor = .x[1],
      Response = .x[2]
    )
  )

stopifnot(
  nrow(model_grid) ==
    choose(
      length(compartment_order),
      2
    )
)

## ------------------------------------------------------------
## 17. Fit all models and adjust for multiple testing
## ------------------------------------------------------------

distance_model_results <- map2_dfr(
  model_grid$Response,
  model_grid$Predictor,

  ~ fit_distance_lmm(
    response_name = .x,
    predictor_name = .y,
    data = bray_pair_table
  )
) %>%
  mutate(
    P_Adjusted_BH =
      p.adjust(
        P_Value,
        method = "BH"
      ),

    Same_Site_P_Adjusted_BH =
      p.adjust(
        Same_Site_P_Value,
        method = "BH"
      ),

    Significant_Raw =
      !is.na(P_Value) &
      P_Value < 0.05,

    Significant_BH =
      !is.na(P_Adjusted_BH) &
      P_Adjusted_BH < 0.05,

    Same_Site_Significant_Raw =
      !is.na(Same_Site_P_Value) &
      Same_Site_P_Value < 0.05,

    Same_Site_Significant_BH =
      !is.na(Same_Site_P_Adjusted_BH) &
      Same_Site_P_Adjusted_BH < 0.05,

    Comparison = paste(
      Response,
      "~",
      Predictor
    )
  ) %>%
  arrange(
    P_Adjusted_BH,
    P_Value
  )

write_csv(
  distance_model_results,
  file.path(
    output_directory,
    "all_pair_dyadic_lmm_results.csv"
  )
)

print(
  distance_model_results %>%
    select(
      Response,
      Predictor,
      Estimate,
      CI_Lower,
      CI_Upper,
      P_Value,
      P_Adjusted_BH,
      Same_Site_Estimate,
      Same_Site_CI_Lower,
      Same_Site_CI_Upper,
      Same_Site_P_Value,
      Same_Site_P_Adjusted_BH,
      N_Pairs,
      N_Within_Site_Pairs,
      N_Between_Site_Pairs,
      N_Sites,
      N_Plants,
      Singular,
      Convergence_Message,
      Model_Status
    ),
  n = Inf
)

## ------------------------------------------------------------
## 18. Human-readable labels and biological comparison groups
## ------------------------------------------------------------

compartment_labels <- c(
  Soil_1m = "Bulk soil (>1 m)",
  Soil_p = "Rhizosphere soil",
  Root = "Root",
  Leaf_EN = "Leaf endophyte",
  Leaf_EP = "Leaf epiphyte",
  Fruit_EN = "Fruit endophyte",
  Fruit_EP = "Fruit epiphyte"
)

compartment_region <- c(
  Soil_1m = "Belowground",
  Soil_p = "Belowground",
  Root = "Belowground",
  Leaf_EN = "Leaf",
  Leaf_EP = "Leaf",
  Fruit_EN = "Fruit",
  Fruit_EP = "Fruit"
)

classify_comparison <- function(
    response,
    predictor
) {
  response_region <-
    unname(
      compartment_region[response]
    )

  predictor_region <-
    unname(
      compartment_region[predictor]
    )

  case_when(
    response_region == "Belowground" &
      predictor_region == "Belowground" ~
      "Belowground pathway",

    response %in% c(
      "Leaf_EN",
      "Leaf_EP"
    ) &
      predictor %in% c(
        "Leaf_EN",
        "Leaf_EP"
      ) ~
      "Leaf pathway",

    response %in% c(
      "Fruit_EN",
      "Fruit_EP"
    ) &
      predictor %in% c(
        "Fruit_EN",
        "Fruit_EP"
      ) ~
      "Fruit pathway",

    response_region %in% c(
      "Leaf",
      "Fruit"
    ) &
      predictor_region %in% c(
        "Leaf",
        "Fruit"
      ) &
      response_region != predictor_region ~
      "Cross-organ: leaf-fruit",

    response_region %in% c(
      "Belowground",
      "Leaf"
    ) &
      predictor_region %in% c(
        "Belowground",
        "Leaf"
      ) &
      response_region != predictor_region ~
      "Cross-organ: belowground-leaf",

    response_region %in% c(
      "Belowground",
      "Fruit"
    ) &
      predictor_region %in% c(
        "Belowground",
        "Fruit"
      ) &
      response_region != predictor_region ~
      "Cross-organ: belowground-fruit",

    TRUE ~ "Other"
  )
}

comparison_category_order <- c(
  "Belowground pathway",
  "Leaf pathway",
  "Fruit pathway",
  "Cross-organ: belowground-leaf",
  "Cross-organ: belowground-fruit",
  "Cross-organ: leaf-fruit"
)

distance_model_results <-
  distance_model_results %>%
  mutate(
    Comparison_Category =
      classify_comparison(
        response = Response,
        predictor = Predictor
      ),

    Comparison_Category =
      factor(
        Comparison_Category,
        levels =
          comparison_category_order
      )
  )

write_csv(
  distance_model_results,
  file.path(
    output_directory,
    "all_pair_dyadic_lmm_results.csv"
  )
)

print(
  distance_model_results %>%
    count(
      Comparison_Category,
      .drop = FALSE
    )
)

if (
  any(
    is.na(
      distance_model_results$
        Comparison_Category
    )
  )
) {
  warning(
    paste(
      "At least one comparison was not",
      "assigned to a category."
    )
  )
}

## ============================================================
## FIGURE 1: COLOR-CODED COEFFICIENT MATRIX
## ============================================================

## ------------------------------------------------------------
## 19. Prepare lower-triangle heatmap data
## ------------------------------------------------------------

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

print(
  coefficient_heatmap
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

## ============================================================
## FIGURE 2: EFFECT-SIZE PLOT BY BIOLOGICAL CATEGORY
## ============================================================

## ------------------------------------------------------------
## 20. Prepare categorized coefficient-plot data
## ------------------------------------------------------------

coefficient_plot_data <-
  distance_model_results %>%
  filter(
    Model_Status == "Success",
    is.finite(Estimate),
    is.finite(CI_Lower),
    is.finite(CI_Upper)
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
      Significant_BH,
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
    Estimate
  ) %>%
  mutate(
    Comparison_Label =
      forcats::fct_inorder(
        Comparison_Label
      )
  )

coefficient_plot <- ggplot(
  coefficient_plot_data,
  aes(
    x = Estimate,
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
      xmin = CI_Lower,
      xmax = CI_Upper
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
        "Associations between",
        "microbiome dissimilarities"
      ),

    subtitle =
      paste(
        "All plant pairs;",
        "models control whether",
        "the pair came from the same site"
      ),

    x =
      "Standardized coefficient (95% CI)",

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

print(
  coefficient_plot
)

ggsave(
  filename = file.path(
    output_directory,
    paste0(
      "all_pair_coefficient_plot_",
      "95CI_by_category.png"
    )
  ),
  plot = coefficient_plot,
  width = 12,
  height = 11,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_directory,
    paste0(
      "all_pair_coefficient_plot_",
      "95CI_by_category.pdf"
    )
  ),
  plot = coefficient_plot,
  width = 12,
  height = 11,
  units = "in"
)

## ============================================================
## FIGURE 3: SAME-SITE EFFECTS
## ============================================================

## ------------------------------------------------------------
## 21. Plot the Same_Site coefficient from every model
## ------------------------------------------------------------

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

print(
  same_site_plot
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

## ------------------------------------------------------------
## 22. Model-diagnostic summary
## ------------------------------------------------------------

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

print(
  model_diagnostics
)

write_csv(
  model_diagnostics,
  file.path(
    output_directory,
    "model_diagnostic_summary.csv"
  )
)

## ============================================================
## SITE EFFECT ON EACH MICROBIOME COMPARTMENT
## ============================================================
##
## Question:
## What proportion of variation in microbial community composition
## is associated with site within each sampled compartment?
##
## For each compartment, this section:
## 1. Calculates Bray-Curtis dissimilarity among samples.
## 2. Runs PERMANOVA:
##
##      Bray-Curtis dissimilarity ~ Site
##
## 3. Extracts the PERMANOVA R2 for Site.
## 4. Tests homogeneity of multivariate dispersion among sites.
## 5. Applies BH correction separately to:
##      - the seven PERMANOVA tests
##      - the seven dispersion tests
## 6. Produces a figure showing site R2 by compartment.
##
## Interpretation:
## PERMANOVA R2 is the fraction of variation in the Bray-Curtis
## distance matrix associated with site.
##
## A significant dispersion test means that sites differ in their
## within-site variability, which should be considered when interpreting
## a significant PERMANOVA result.
## ============================================================

## ------------------------------------------------------------
## 23. Function to run site PERMANOVA for one compartment
## ------------------------------------------------------------

run_site_permanova <- function(
    compartment_name,
    metadata_table,
    abundance_matrix,
    permutations = 9999
) {

  compartment_metadata <- metadata_table %>%
    filter(
      Compartment == compartment_name
    ) %>%
    arrange(
      Site,
      suppressWarnings(
        as.numeric(Plant)
      )
    )

  compartment_abundance <- abundance_matrix[
    compartment_metadata$Sample_Index,
    ,
    drop = FALSE
  ]

  # Remove samples with no reads.
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

  # Remove taxa absent from all retained samples.
  if (nrow(compartment_abundance) > 0) {
    keep_taxa <- colSums(
      compartment_abundance,
      na.rm = TRUE
    ) > 0

    compartment_abundance <- compartment_abundance[
      ,
      keep_taxa,
      drop = FALSE
    ]
  }

  number_of_samples <- nrow(
    compartment_metadata
  )

  number_of_sites <- dplyr::n_distinct(
    compartment_metadata$Site
  )

  empty_result <- function(
      status_message
  ) {
    tibble(
      Compartment = compartment_name,
      N_Samples = number_of_samples,
      N_Sites = number_of_sites,

      Site_DF = NA_real_,
      Residual_DF = NA_real_,

      Site_SumOfSqs = NA_real_,
      Residual_SumOfSqs = NA_real_,

      Site_R2 = NA_real_,
      Residual_R2 = NA_real_,

      PERMANOVA_F = NA_real_,
      PERMANOVA_P = NA_real_,

      Dispersion_F = NA_real_,
      Dispersion_P = NA_real_,

      Analysis_Status = status_message
    )
  }

  if (
    number_of_samples < 3 ||
    number_of_sites < 2 ||
    ncol(compartment_abundance) < 1
  ) {
    return(
      empty_result(
        "Insufficient samples, sites, or taxa"
      )
    )
  }

  compartment_metadata <- compartment_metadata %>%
    mutate(
      Site = droplevels(
        factor(Site)
      )
    )

  site_sample_counts <- table(
    compartment_metadata$Site
  )

  if (any(site_sample_counts < 2)) {
    warning(
      paste(
        "At least one site contains fewer than two samples in",
        compartment_name,
        "- dispersion results may be unreliable."
      )
    )
  }

  # Calculate Bray-Curtis distances.
  bray_distance <- tryCatch(
    vegan::vegdist(
      compartment_abundance,
      method = "bray"
    ),
    error = function(e) e
  )

  if (inherits(bray_distance, "error")) {
    return(
      empty_result(
        paste(
          "Bray-Curtis error:",
          bray_distance$message
        )
      )
    )
  }

  if (
    any(
      !is.finite(
        as.vector(bray_distance)
      )
    )
  ) {
    return(
      empty_result(
        "Non-finite Bray-Curtis distances"
      )
    )
  }

  # Minimal metadata frame aligned with the distance matrix.
  model_data <- data.frame(
    Site = compartment_metadata$Site
  )

  set.seed(12345)

  permanova_result <- tryCatch(
    vegan::adonis2(
      bray_distance ~ Site,
      data = model_data,
      permutations = permutations,
      by = "terms"
    ),
    error = function(e) e
  )

  if (inherits(permanova_result, "error")) {
    return(
      empty_result(
        paste(
          "PERMANOVA error:",
          permanova_result$message
        )
      )
    )
  }

  permanova_table <- as.data.frame(
    permanova_result
  )

  permanova_table$Term <- rownames(
    permanova_table
  )

  # Depending on vegan version/output settings, the tested term may be
  # labeled either "Site" or "Model".
  site_row <- permanova_table[
    permanova_table$Term %in% c(
      "Site",
      "Model"
    ),
    ,
    drop = FALSE
  ]

  residual_row <- permanova_table[
    permanova_table$Term == "Residual",
    ,
    drop = FALSE
  ]

  if (
    nrow(site_row) != 1 ||
    nrow(residual_row) != 1
  ) {
    return(
      empty_result(
        paste0(
          "Could not identify Site/Model and Residual rows. Terms returned: ",
          paste(
            permanova_table$Term,
            collapse = ", "
          )
        )
      )
    )
  }

  permanova_p_column <- grep(
    "^Pr\\(",
    names(permanova_table),
    value = TRUE
  )

  if (length(permanova_p_column) != 1) {
    return(
      empty_result(
        paste0(
          "Could not identify PERMANOVA p-value column. Columns returned: ",
          paste(
            names(permanova_table),
            collapse = ", "
          )
        )
      )
    )
  }

  # Test equality of multivariate dispersion among sites.
  dispersion_model <- tryCatch(
    vegan::betadisper(
      bray_distance,
      group = model_data$Site,
      type = "centroid",
      bias.adjust = TRUE
    ),
    error = function(e) e
  )

  if (inherits(dispersion_model, "error")) {
    return(
      empty_result(
        paste(
          "Dispersion-model error:",
          dispersion_model$message
        )
      )
    )
  }

  set.seed(12345)

  dispersion_result <- tryCatch(
    vegan::permutest(
      dispersion_model,
      permutations = permutations
    ),
    error = function(e) e
  )

  if (inherits(dispersion_result, "error")) {
    return(
      empty_result(
        paste(
          "Dispersion-test error:",
          dispersion_result$message
        )
      )
    )
  }

  dispersion_table <- as.data.frame(
    dispersion_result$tab
  )

  dispersion_p_column <- grep(
    "^Pr\\(",
    names(dispersion_table),
    value = TRUE
  )

  if (length(dispersion_p_column) != 1) {
    return(
      empty_result(
        paste0(
          "Could not identify dispersion p-value column. Columns returned: ",
          paste(
            names(dispersion_table),
            collapse = ", "
          )
        )
      )
    )
  }

  tibble(
    Compartment = compartment_name,
    N_Samples = number_of_samples,
    N_Sites = number_of_sites,

    Site_DF = as.numeric(
      site_row$Df
    ),

    Residual_DF = as.numeric(
      residual_row$Df
    ),

    Site_SumOfSqs = as.numeric(
      site_row$SumOfSqs
    ),

    Residual_SumOfSqs = as.numeric(
      residual_row$SumOfSqs
    ),

    Site_R2 = as.numeric(
      site_row$R2
    ),

    Residual_R2 = as.numeric(
      residual_row$R2
    ),

    PERMANOVA_F = as.numeric(
      site_row$F
    ),

    PERMANOVA_P = as.numeric(
      site_row[[permanova_p_column]]
    ),

    Dispersion_F = as.numeric(
      dispersion_table$F[1]
    ),

    Dispersion_P = as.numeric(dispersion_table[[dispersion_p_column]][1]),

    Analysis_Status = "Success"
  )
}

## ------------------------------------------------------------
## 24. Run site analysis for all seven compartments
## ------------------------------------------------------------

site_permanova_results <- purrr::map_dfr(
  compartment_order,
  run_site_permanova,
  metadata_table = metadata_matched,
  abundance_matrix = taxa_abundance_matched,
  permutations = 9999
) %>%
  mutate(
    PERMANOVA_P_Adjusted_BH = p.adjust(
      PERMANOVA_P,
      method = "BH"
    ),

    Dispersion_P_Adjusted_BH = p.adjust(
      Dispersion_P,
      method = "BH"
    ),

    PERMANOVA_Significant_BH =
      !is.na(PERMANOVA_P_Adjusted_BH) &
      PERMANOVA_P_Adjusted_BH < 0.05,

    Dispersion_Significant_BH =
      !is.na(Dispersion_P_Adjusted_BH) &
      Dispersion_P_Adjusted_BH < 0.05
  )

print(
  site_permanova_results %>%
    select(
      Compartment,
      N_Samples,
      N_Sites,
      Site_R2,
      PERMANOVA_P,
      PERMANOVA_P_Adjusted_BH,
      Dispersion_P,
      Dispersion_P_Adjusted_BH,
      Analysis_Status
    ),
  n = Inf,
  width = Inf
)

## ------------------------------------------------------------
## 25. Add readable compartment labels
## ------------------------------------------------------------

# compartment_labels already exists earlier in the script. ###
site_permanova_results <- site_permanova_results %>%
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

    PERMANOVA_Result = if_else(
      PERMANOVA_Significant_BH,
      "Site effect: BH-adjusted p < 0.05",
      "Site effect: BH-adjusted p >= 0.05"
    ),

    Dispersion_Result = if_else(
      Dispersion_Significant_BH,
      "Dispersion differs among sites",
      "No significant dispersion difference"
    )
  )

## ------------------------------------------------------------
## 26. Save site PERMANOVA results
## ------------------------------------------------------------

write_csv(
  site_permanova_results,
  file.path(
    output_directory,
    "site_PERMANOVA_R2_by_compartment.csv"
  )
)

print(
  site_permanova_results %>%
    select(
      Compartment_Label,
      N_Samples,
      N_Sites,
      Site_R2,
      PERMANOVA_F,
      PERMANOVA_P,
      PERMANOVA_P_Adjusted_BH,
      Dispersion_F,
      Dispersion_P,
      Dispersion_P_Adjusted_BH,
      Analysis_Status
    ),
  n = Inf,
  width = Inf
)

## ============================================================
## FIGURE 4: SITE R2 BY COMPARTMENT
## ============================================================

## ------------------------------------------------------------
## 27. Prepare site-R2 plot data
## ------------------------------------------------------------

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

site_r2_upper_limit <- max(
  0.10,
  max(
    site_r2_plot_data$Site_R2,
    na.rm = TRUE
  ) * 1.30
)

## ------------------------------------------------------------
## 28. Plot site PERMANOVA R2
## ------------------------------------------------------------

site_r2_plot <- ggplot(
  site_r2_plot_data,
  aes(
    x = Site_R2,
    y = Compartment_Label,
    color = PERMANOVA_Result,
    shape = Dispersion_Result
  )
) +
  geom_segment(
    aes(
      x = 0,
      xend = Site_R2,
      y = Compartment_Label,
      yend = Compartment_Label
    ),
    color = "grey75",
    linewidth = 0.8
  ) +
  geom_point(
    size = 4
  ) +
  geom_text(
    aes(
      label = sprintf(
        "R² = %.2f",
        Site_R2
      )
    ),
    hjust = -0.15,
    color = "black",
    size = 3.7
  ) +
  scale_color_manual(
    values = c(
      "Site effect: BH-adjusted p < 0.05" =
        "#B2182B",

      "Site effect: BH-adjusted p >= 0.05" =
        "#4D4D4D"
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
    limits = c(
      0,
      site_r2_upper_limit
    ),
    expand = expansion(
      mult = c(
        0,
        0.02
      )
    )
  ) +
  labs(
    title =
      "Strength of site effects across microbiome compartments",

    subtitle = paste(
      "Points show the proportion of Bray-Curtis variation",
      "associated with site; significance uses BH-adjusted",
      "PERMANOVA p-values"
    ),

    x = expression(
      "Site PERMANOVA " * R^2
    ),

    y = NULL,

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
        size = 11
      ),

    plot.title =
      element_text(
        face = "bold"
      ),

    legend.position =
      "top"
  )

print(
  site_r2_plot
)

ggsave(
  filename = file.path(
    output_directory,
    "site_PERMANOVA_R2_by_compartment.png"
  ),
  plot = site_r2_plot,
  width = 10,
  height = 6,
  units = "in",
  dpi = 400,
  bg = "white"
)

ggsave(
  filename = file.path(
    output_directory,
    "site_PERMANOVA_R2_by_compartment.pdf"
  ),
  plot = site_r2_plot,
  width = 10,
  height = 6,
  units = "in"
)


## ============================================================
## FIGURE 5: WITHIN-SITE VS AMONG-SITE BRAY-CURTIS DISTRIBUTIONS
## ============================================================

library(ggdist)

## ------------------------------------------------------------
## 29. Prepare staggered half-violin data
## ------------------------------------------------------------

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

## ------------------------------------------------------------
## 30. Draw staggered upper-half violins
## ------------------------------------------------------------

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

print(
  within_among_plot
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

## ============================================================
## SAMPLE-LEVEL READ DEPTH AND RICHNESS CHECKS
## ============================================================

## ------------------------------------------------------------
## 31. Calculate read depth and observed genus richness
## ------------------------------------------------------------

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

## ------------------------------------------------------------
## 32. Richness plot
## ------------------------------------------------------------

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

print(
  richness_plot
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

## ------------------------------------------------------------
## 33. Read-depth plot
## ------------------------------------------------------------

read_depth_plot <- ggplot(
  sample_richness,
  aes(
    x = Total_Reads,
    y = Compartment_Label
  )
) +
  geom_boxplot(
    width = 0.6,
    outlier.alpha = 0.7
  ) +
  scale_x_log10(
    labels = scales::label_comma()
  ) +
  labs(
    title = "Sequencing depth across microbiome compartments",
    x = "Total reads (log10 scale)",
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

print(
  read_depth_plot
)

ggsave(
  filename = file.path(
    output_directory,
    "total_reads_by_compartment.png"
  ),
  plot = read_depth_plot,
  width = 8,
  height = 5.5,
  units = "in",
  dpi = 400,
  bg = "white"
)

## ------------------------------------------------------------
## 34. Final reproducibility information
## ------------------------------------------------------------

capture.output(
  sessionInfo(),
  file = file.path(
    output_directory,
    "sessionInfo.txt"
  )
)

message(
  "Analysis complete. Results and figures were saved to: ",
  normalizePath(
    output_directory,
    mustWork = FALSE
  )
)


## ============================================================
## DOES THE WITHIN-VS-AMONG SITE GAP DIFFER BY COMPARTMENT?
## ============================================================
##
## Question:
## Does the increase in Bray-Curtis dissimilarity from within-site
## pairs to among-site pairs differ among microbiome compartments?
##
## The key test is:
##
##     Compartment × Among_Site
##
## Among_Site is coded:
##     0 = plants from the same site
##     1 = plants from different sites
##
## Therefore, the Among_Site slope within each compartment is:
##
##     mean among-site dissimilarity - mean within-site dissimilarity
##
## A small slope means that within-site pairs are nearly as dissimilar
## as among-site pairs.
## ============================================================

library(emmeans)

## ------------------------------------------------------------
## 1. Prepare the dyadic repeated-measures data
## ------------------------------------------------------------

site_gap_data <- bray_long %>%
  filter(
    is.finite(Bray_Curtis)
  ) %>%
  mutate(
    Compartment = factor(
      Compartment,
      levels = compartment_order
    ),
    
    # FALSE = within site; TRUE = among sites
    Among_Site = as.numeric(
      !Same_Site
    ),
    
    Site_1 = factor(Site_1),
    Site_2 = factor(Site_2),
    
    Plant_1_ID = factor(Plant_1_ID),
    Plant_2_ID = factor(Plant_2_ID),
    
    # Each plant pair is repeatedly measured across compartments.
    Pair_ID = factor(
      paste(
        Plant_1_ID,
        Plant_2_ID,
        sep = "__"
      )
    )
  )

# Confirm that complete plant pairs are represented across compartments. ###
site_gap_data %>%
  count(
    Pair_ID,
    name = "N_Compartments"
  ) %>%
  count(
    N_Compartments
  ) %>%
  print(n = Inf)

## ------------------------------------------------------------
## 2. Fit model without compartment-specific site gaps
## ------------------------------------------------------------

site_gap_null_model <- lmer(
  Bray_Curtis ~
    Compartment +
    Among_Site +
    (1 | Site_1) +
    (1 | Site_2) +
    (1 | Plant_1_ID) +
    (1 | Plant_2_ID) +
    (1 | Pair_ID),
  
  data = site_gap_data,
  REML = FALSE,
  
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(
      maxfun = 200000
    )
  )
)

## ------------------------------------------------------------
## 3. Fit interaction model
## ------------------------------------------------------------

site_gap_interaction_model <- lmer(
  Bray_Curtis ~
    Compartment * Among_Site +
    (1 | Site_1) +
    (1 | Site_2) +
    (1 | Plant_1_ID) +
    (1 | Plant_2_ID) +
    (1 | Pair_ID),
  
  data = site_gap_data,
  REML = FALSE,
  
  control = lmerControl(
    optimizer = "bobyqa",
    optCtrl = list(
      maxfun = 200000
    )
  )
)

print(
  summary(
    site_gap_interaction_model
  )
)

## ------------------------------------------------------------
## 4. Omnibus interaction test
## ------------------------------------------------------------
##
## This likelihood-ratio test asks whether allowing a different
## among-vs-within gap for each compartment improves model fit.

site_gap_interaction_test <- anova(
  site_gap_null_model,
  site_gap_interaction_model
)

print(
  site_gap_interaction_test
)

write.csv(
  as.data.frame(
    site_gap_interaction_test
  ),
  file.path(
    output_directory,
    "compartment_by_site_gap_interaction_test.csv"
  ),
  row.names = TRUE
)

## ------------------------------------------------------------
## 5. Estimate among-minus-within gap for each compartment
## ------------------------------------------------------------
##
## Because Among_Site is numeric, emtrends extracts its slope separately
## for each compartment. A positive slope means among-site pairs are more
## dissimilar than within-site pairs.

site_gap_estimates <- emmeans::emtrends(
  site_gap_interaction_model,
  specs = ~ Compartment,
  var = "Among_Site",
  lmer.df = "asymptotic"
)

site_gap_results <- as.data.frame(
  summary(
    site_gap_estimates,
    infer = c(TRUE, TRUE),
    level = 0.95
  )
) %>%
  as_tibble() %>%
  rename(
    Among_Minus_Within = Among_Site.trend,
    Standard_Error = SE,
    Degrees_of_Freedom = df,
    CI_Lower = asymp.LCL,
    CI_Upper = asymp.UCL,
    P_Value = p.value
  ) %>%
  mutate(
    P_Adjusted_BH = p.adjust(
      P_Value,
      method = "BH"
    ),
    
    Significant_BH =
      !is.na(P_Adjusted_BH) &
      P_Adjusted_BH < 0.05,
    
    Compartment_Label = unname(
      compartment_labels[
        as.character(Compartment)
      ]
    )
  )

print(
  site_gap_results,
  n = Inf,
  width = Inf
)

write_csv(
  site_gap_results,
  file.path(
    output_directory,
    "among_minus_within_site_gap_by_compartment.csv"
  )
)

## ------------------------------------------------------------
## 6. Compare site gaps among compartments
## ------------------------------------------------------------
##
## These pairwise contrasts ask whether, for example, the site gap in
## fruit endophytes is smaller than the site gap in leaf endophytes.

site_gap_pairwise <- pairs(
  site_gap_estimates,
  adjust = "BH"
)

site_gap_pairwise_results <- as.data.frame(
  summary(
    site_gap_pairwise,
    infer = c(TRUE, TRUE),
    level = 0.95
  )
) %>%
  as_tibble()

print(
  site_gap_pairwise_results,
  n = Inf,
  width = Inf
)

write_csv(
  site_gap_pairwise_results,
  file.path(
    output_directory,
    "pairwise_comparisons_of_site_gaps.csv"
  )
)

## ============================================================
## FIGURE: AMONG-MINUS-WITHIN SITE GAP
## ============================================================

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

print(
  site_gap_plot
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