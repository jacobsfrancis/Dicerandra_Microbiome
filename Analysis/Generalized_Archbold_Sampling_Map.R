# Map of Archbold sampling locations with Florida inset ####
#
# This map uses representative site coordinates calculated from the median
# latitude and longitude of sampled plants within each site.
#
# Five sites are highlighted as included in the present study.
# Three additional sampled sites are shown in charcoal.
#
# The Florida inset identifies the regional location of Archbold Biological
# Station, while the main panel shows the local station boundary.

library(sf)
library(tidyverse)
library(maps)
library(ggspatial)
library(patchwork)

# 0) File paths ####

boundary_zip <- "Archbold_Biological_Station_-_Current_Boundaries.zip"

metadata_file <- "../Data/FPCP_metadata.csv"

output_directory <- file.path(
  "../Output/",
  "Poster_Common_5_Sites"
)

dir.create(
  output_directory,
  showWarnings = FALSE,
  recursive = TRUE
)

# 1) Extract and read the ABS boundary ####

boundary_directory <- file.path(
  tempdir(),
  "archbold_boundary"
)

dir.create(
  boundary_directory,
  showWarnings = FALSE,
  recursive = TRUE
)

unzip(
  zipfile = boundary_zip,
  exdir = boundary_directory
)

boundary_shapefile <- list.files(
  boundary_directory,
  pattern = "\\.shp$",
  full.names = TRUE
)

if (length(boundary_shapefile) != 1) {
  stop(
    "Expected exactly one shapefile in the boundary archive."
  )
}

archbold_boundary <- st_read(
  boundary_shapefile,
  quiet = TRUE
) %>%
  st_make_valid() %>%
  st_transform(
    3857
  )

# Combine all boundary polygons into one region.

archbold_union <- archbold_boundary %>%
  summarise(
    geometry = st_union(geometry)
  )

# Create a longitude-latitude copy for the Florida inset.

archbold_lonlat <- archbold_union %>%
  st_transform(
    4326
  )

archbold_centroid_lonlat <- archbold_lonlat %>%
  st_point_on_surface()

archbold_coordinates <- st_coordinates(
  archbold_centroid_lonlat
)

archbold_longitude <- archbold_coordinates[1, "X"]
archbold_latitude <- archbold_coordinates[1, "Y"]

# 2) Read sampling coordinates from metadata ####

metadata <- read_csv(
  metadata_file,
  show_col_types = FALSE
)

included_sites <- c(
  1,
  2,
  4,
  5,
  8
)

# Samples from the same plant occur repeatedly because multiple microbial
# compartments were collected. Retain one coordinate record per plant before
# calculating a representative coordinate for each site.

site_coordinates <- metadata %>%
  transmute(
    Site = as.integer(Sample_Site),
    Plant = Sample_Plant,
    Latitude = as.numeric(Sample_Lat),
    Longitude = as.numeric(Sample_Long)
  ) %>%
  filter(
    !is.na(Site),
    !is.na(Plant),
    !is.na(Latitude),
    !is.na(Longitude)
  ) %>%
  distinct(
    Site,
    Plant,
    Latitude,
    Longitude
  ) %>%
  group_by(
    Site
  ) %>%
  summarise(
    Longitude = median(
      Longitude,
      na.rm = TRUE
    ),
    Latitude = median(
      Latitude,
      na.rm = TRUE
    ),
    N_Plants = n_distinct(Plant),
    .groups = "drop"
  ) %>%
  mutate(
    Site_Status = if_else(
      Site %in% included_sites,
      "Included in this study",
      "Additional population"
    ),
    Site_Status = factor(
      Site_Status,
      levels = c(
        "Included in this study",
        "Additional population"
      )
    )
  )

# Confirm that eight sites were recovered.

if (nrow(site_coordinates) != 8) {
  stop(
    paste(
      "Expected coordinates for eight sites, but found",
      nrow(site_coordinates)
    )
  )
}

# Convert longitude and latitude to an sf object.

sampling_locations_lonlat <- st_as_sf(
  site_coordinates,
  coords = c(
    "Longitude",
    "Latitude"
  ),
  crs = 4326,
  remove = FALSE
)

# Transform sampling points into the boundary coordinate system.

sampling_locations <- sampling_locations_lonlat %>%
  st_transform(
    st_crs(archbold_union)
  )

# Check whether representative site points fall within the ABS boundary.

inside_archbold <- lengths(
  st_within(
    sampling_locations,
    archbold_union
  )
) > 0

coordinate_check <- site_coordinates %>%
  mutate(
    Inside_Archbold = inside_archbold
  )

print(coordinate_check)

if (any(!inside_archbold)) {
  warning(
    "One or more site coordinates fall outside the mapped Archbold boundary."
  )
}

# 3) Define site-status colors ####

site_status_colors <- c(
  "Included in this study" = "#995C8F",
  "Additional population" = "#4A4A4A"
)

# 4) Prepare the Florida inset ####

florida_map <- ggplot2::map_data(
  "state"
) %>%
  filter(
    region %in% c(
      "florida",
      "georgia",
      "alabama"
    )
  )

# The rectangle marks the approximate regional position of Archbold rather
# than the precise station footprint.

inset_box_half_width_degrees <- 0.35
inset_box_half_height_degrees <- 0.30

archbold_inset_box <- tibble(
  xmin = archbold_longitude - inset_box_half_width_degrees,
  xmax = archbold_longitude + inset_box_half_width_degrees,
  ymin = archbold_latitude - inset_box_half_height_degrees,
  ymax = archbold_latitude + inset_box_half_height_degrees
)

# 5) Define the map theme ####

map_theme <- theme_void(
  base_size = 18,
  base_family = "Helvetica"
) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0,
      margin = margin(
        b = 3
      )
    ),
    plot.subtitle = element_text(
      hjust = 0,
      margin = margin(
        b = 5
      )
    ),
    legend.title = element_blank(),
    legend.text = element_text(
      size = 16
    ),
    legend.position = "bottom",
    legend.direction = "vertical",
    legend.justification = "left",
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

# 6) Draw the Florida inset ####

florida_inset_plot <- ggplot() +
  
  geom_polygon(
    data = florida_map,
    aes(
      x = long,
      y = lat,
      group = group
    ),
    fill = "#F2F2F2",
    color = "#4A4A4A",
    linewidth = 0.45
  ) +
  
  geom_rect(
    data = archbold_inset_box,
    aes(
      xmin = xmin,
      xmax = xmax,
      ymin = ymin,
      ymax = ymax
    ),
    inherit.aes = FALSE,
    fill = NA,
    color = "#995C8F",
    linewidth = 0.9,
    linetype = "dashed"
  ) +
  
  geom_point(
    aes(
      x = archbold_longitude,
      y = archbold_latitude
    ),
    shape = 21,
    size = 4.5,
    stroke = 0.9,
    fill = "#995C8F",
    color = "#262626"
  ) +
  
  annotate(
    geom = "text",
    x = archbold_longitude - 0.10,
    y = archbold_latitude + 0.26,
    label = "Archbold Biological Station",
    hjust = 1,
    vjust = 0,
    family = "Helvetica",
    fontface = "bold",
    size = 14 / .pt,
    color = "#262626"
  ) +
  
  coord_quickmap(
    xlim = c(
      -87.8,
      -79.6
    ),
    ylim = c(
      24.3,
      31.1
    ),
    expand = FALSE
  ) +
  
  labs(
    title = "Florida"
  ) +
  
  map_theme +
  
  theme(
    legend.position = "none",
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    panel.border = element_rect(
      fill = NA,
      color = "#4A4A4A",
      linewidth = 0.5
    )
  )

# 7) Draw the Archbold map ####

archbold_map <- ggplot() +
  
  geom_sf(
    data = archbold_union,
    fill = "#F2F2F2",
    color = "#262626",
    linewidth = 0.9
  ) +
  
  geom_sf(
    data = sampling_locations,
    aes(
      fill = Site_Status
    ),
    shape = 21,
    size = 5.8,
    stroke = 1,
    color = "#262626"
  ) +
  
  scale_fill_manual(
    values = site_status_colors,
    labels = c(
      "Included in this study" =
        "Included in this study (n = 5)",
      "Additional population" =
        "Not included in this study (n = 3)"
    ),
    drop = FALSE
  ) +
  
  coord_sf(
    datum = NA,
    expand = FALSE
  ) +
  
  labs(
    title = "Archbold Biological Station",
    subtitle = "Points show representative coordinates for the eight sampled sites."
  ) +
  
  annotation_scale(
    location = "bl",
    width_hint = 0.25,
    text_cex = 0.9,
    line_width = 0.8,
    pad_x = unit(
      0.25,
      "cm"
    ),
    pad_y = unit(
      0.25,
      "cm"
    )
  ) +
  
  map_theme

# 8) Assemble the inset and Archbold map ####

sampling_location_map <-
  florida_inset_plot +
  archbold_map +
  patchwork::plot_layout(
    widths = c(
      0.72,
      1.55
    )
  )

sampling_location_map

# 9) Save the map ####

ggsave(
  filename = file.path(
    output_directory,
    "archbold_sampling_locations_with_florida_inset.png"
  ),
  plot = sampling_location_map,
  width = 6,
  height = 5,
  units = "in",
  dpi = 600,
  bg = "transparent"
)

ggsave(
  filename = file.path(
    output_directory,
    "archbold_sampling_locations_with_florida_inset.pdf"
  ),
  plot = sampling_location_map,
  width = 6,
  height = 5,
  units = "in",
  bg = "transparent"
)

# 10) Save summarized site coordinates ####

write_csv(
  site_coordinates,
  file.path(
    output_directory,
    "archbold_site_coordinates_summary.csv"
  )
)