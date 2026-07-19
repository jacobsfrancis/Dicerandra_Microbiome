# 2) Run statistical analyses ####

library(tidyverse)
library(vegan)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(emmeans)

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

# 2) Test coordination among microbiome compartments ####

## 2a) Fit an example coordination model ####

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

## 2b) Define the coordination-model function ####

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

## 2c) Define all unique compartment comparisons ####
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

## 2d) Fit coordination models and control false discovery rate ####

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

## 2e) Label and classify compartment comparisons ####

compartment_labels <- c(
  Soil_1m = "Bulk soil (>1 m)",
  Soil_p = "Rhizosphere soil",
  Root = "Root epiphyte",
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
      "Leaf-fruit",

    response_region %in% c(
      "Belowground",
      "Leaf"
    ) &
      predictor_region %in% c(
        "Belowground",
        "Leaf"
      ) &
      response_region != predictor_region ~
      "Belowground-leaf",

    response_region %in% c(
      "Belowground",
      "Fruit"
    ) &
      predictor_region %in% c(
        "Belowground",
        "Fruit"
      ) &
      response_region != predictor_region ~
      "Belowground-fruit",

    TRUE ~ "Other"
  )
}

comparison_category_order <- c(
  "Belowground pathway",
  "Leaf pathway",
  "Fruit pathway",
  "Belowground-leaf",
  "Belowground-fruit",
  "Leaf-fruit"
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

# 3) Quantify site effects within compartments ####
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

## 3a) Define the site PERMANOVA and dispersion tests ####

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

## 3b) Run site analyses across compartments ####

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

## 3c) Label site-analysis results ####

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

## 3d) Save site-analysis results ####

write_csv(
  site_permanova_results,
  file.path(
    output_directory,
    "site_PERMANOVA_R2_by_compartment.csv"
  )
)

# 4) Test whether site-associated similarity differs among compartments ####
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

library(emmeans)

## 4a) Prepare repeated dyadic observations ####

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
  ) 
## 4b) Fit the additive model ####

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

## 4c) Fit the compartment-by-site-status interaction ####

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

## 4d) Test the omnibus interaction ####
##
## This likelihood-ratio test asks whether allowing a different
## among-vs-within gap for each compartment improves model fit.

site_gap_interaction_test <- anova(
  site_gap_null_model,
  site_gap_interaction_model
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

## 4e) Estimate compartment-specific site gaps ####
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

write_csv(
  site_gap_results,
  file.path(
    output_directory,
    "among_minus_within_site_gap_by_compartment.csv"
  )
)

## 4f) Compare site gaps among compartments ####
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

write_csv(
  site_gap_pairwise_results,
  file.path(
    output_directory,
    "pairwise_comparisons_of_site_gaps.csv"
  )
)

# 5) Save analysis objects for figure scripts ####

analysis_results <- list(
  distance_model_results = distance_model_results,
  site_permanova_results = site_permanova_results,
  site_gap_results = site_gap_results,
  site_gap_pairwise_results = site_gap_pairwise_results,
  site_gap_interaction_test = site_gap_interaction_test,
  compartment_labels = compartment_labels,
  compartment_region = compartment_region,
  comparison_category_order = comparison_category_order
)

saveRDS(
  analysis_results,
  file.path(intermediate_directory, "analysis_results.rds")
)

capture.output(
  sessionInfo(),
  file = file.path(output_directory, "sessionInfo.txt")
)
