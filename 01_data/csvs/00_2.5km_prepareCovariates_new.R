# =============================================================================
# COVARIATE AGGREGATION TO 2.5KM GRID
# Processes climate, land, and observational covariates and builds final
# standardised covariate datasets for modelling
# =============================================================================

# ---- CONFIG -----------------------------------------------------------------

# --- Project paths
grid_1km_path       <- "01_data/grids/grid_clipped_1km.gpkg"
grid_2point5km_path <- "01_data/grids/grid_clipped_2point5km.gpkg"
scotland_path       <- "01_data/gridExtractions/scotlandBoundary/scotland_boundary.shp"

# --- Raw input data (local/OneDrive)
livestock_rast_path <- "~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/livestockDensity/GLW4-2020.D-DA.GLEAM3-ALL-LU.tif"
elevation_rast_path <- "~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/elevation_DEM/scotland_elevation.tif"
landcover_rast_path <- "01_data/gridExtractions/landcover_27700.tif"
buildings_path      <- "~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/scotlandOSMBuildings_Others.shp/gis_osm_buildings_a_free_1.shp"
poi_path            <- "~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/POI/poi_uk.gpkg"

# --- Processed climate inputs (from temperature/rainfall extraction scripts)
tas_1km_path  <- "01_data/csvs/tasmean_grid_21d_celsius.csv"
rainfall_1km_path <- "01_data/csvs/rainfall_grid_new.csv"
cs_path          <- "01_data/csvs/cs_df.csv"

# --- Output paths
out_dir <- "01_data/covariates"

out_dist_matrix  <- "01_data/covariates/dist_matrix_2point5km_meters.Rdata"
out_temp         <- file.path(out_dir, "mean_temp_21day_2point5km.csv")
out_rainfall     <- file.path(out_dir, "rainfall_28day_2point5km.csv")
out_livestock    <- file.path(out_dir, "livestockDensity_2point5km.csv")
out_elevation    <- file.path(out_dir, "elevation_2point5km.csv")
out_landcover    <- file.path(out_dir, "landcover_2point5km.csv")
out_landcover_grp <- file.path(out_dir, "landcover_grouped_2point5km.csv")
out_buildings    <- file.path(out_dir, "buildings_2point5km.csv")
out_poi          <- file.path(out_dir, "poi_grouped_2point5km.csv")
out_reports      <- file.path(out_dir, "cs_reports_2point5km.csv")
out_z_land       <- file.path(out_dir, "z_land_2point5km.csv")
out_z_climate    <- file.path(out_dir, "z_climate_2point5km.csv")

# --- Parameters
landcover_num_cores <- 10
landcover_chunk_size <- 2000

# Landcover class lookup
landcover_lookup <- c(
  `0`  = "Unnamed",           `1`  = "Deciduous_woodland", `2`  = "Coniferous_woodland",
  `3`  = "Arable",            `4`  = "Improved_grassland", `5`  = "Neutral_grassland",
  `6`  = "Calcareous_grassland", `7` = "Acid_grassland",   `8`  = "Fen",
  `9`  = "Heather",           `10` = "Heather_grassland",  `11` = "Bog",
  `12` = "Inland_rock",       `13` = "Saltwater",          `14` = "Freshwater",
  `15` = "Supralittoral_rock", `16` = "Supralittoral_sediment", `17` = "Littoral_rock",
  `18` = "Littoral_sediment", `19` = "Saltmarsh",          `20` = "Urban",
  `21` = "Suburban"
)

# ---- LIBRARIES --------------------------------------------------------------

library(sf)
library(terra)
library(readr)
library(dplyr)
library(exactextractr)
library(tidyr)
library(future.apply)

# ---- LOAD BASE GRIDS --------------------------------------------------------

message("Loading grids...")
grid_1km       <- st_read(grid_1km_path)
grid_2point5km <- st_read(grid_2point5km_path)
scotland       <- st_read(scotland_path) |> vect()

grid_2point5km$area <- as.numeric(st_area(grid_2point5km)) / 1e6  # km²

# 1km -> 2.5km mapping (used for climate aggregation)
mapping_1km_to_2point5km <- st_join(
  grid_1km       |> dplyr::select(grid_id),
  grid_2point5km |> dplyr::select(grid_id),
  left = FALSE
) |>
  st_drop_geometry() |>
  rename(grid_id_1km = grid_id.x, grid_id_2point5km = grid_id.y)


# ---- DISTANCE MATRIX --------------------------------------------------------

message("Computing distance matrix...")
centroids_2point5km <- st_centroid(grid_2point5km)
dist_matrix_2point5km_meters <- as.matrix(st_distance(centroids_2point5km))
save(dist_matrix_2point5km_meters, file = out_dist_matrix)
message("  Saved: ", out_dist_matrix)


# ---- HELPER: WKT GRID EXPORT ------------------------------------------------

grid_2point5km_wkt <- grid_2point5km |>
  mutate(geometry = st_as_text(geom))
write_csv(grid_2point5km_wkt, "01_data/grids/grid_clipped_2point5km_wkt.csv")
write_csv(grid_2point5km_wkt, "01_data/csvs/grid_clipped_2point5km_wkt.csv")


# ---- 1. TEMPERATURE ---------------------------------------------------------

message("Processing minimum temperature data...")

temp_1km <- read_csv(tas_1km_path)

temp_2point5km <- temp_1km |>
  left_join(mapping_1km_to_2point5km, by = c("grid_id" = "grid_id_1km")) |>
  group_by(grid_id_2point5km, date) |>
  summarise(
    min_temp_21d_celsius = mean(mean_min_temp_21d_celsius, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    min_temp_21d_celsius = as.numeric(min_temp_21d_celsius),
    z_temp_min = scale(min_temp_21d_celsius)[, 1]
  ) |>
  rename(grid_id = grid_id_2point5km)

write_csv(temp_2point5km, out_temp)
message("  Saved: ", out_temp)


# ---- 2. RAINFALL ------------------------------------------------------------

message("Processing rainfall data...")

rainfall_1km <- read_csv(rainfall_1km_path)

rainfall_2point5km <- rainfall_1km |>
  left_join(mapping_1km_to_2point5km, by = c("grid_id" = "grid_id_1km")) |>
  group_by(grid_id_2point5km, date) |>
  summarise(rainfall_28d = mean(rainfall_28d, na.rm = TRUE), .groups = "drop") |>
  mutate(
    rainfall_28d = as.numeric(rainfall_28d),
    z_rain = scale(rainfall_28d)[, 1]
  ) |>
  rename(grid_id = grid_id_2point5km)

write_csv(rainfall_2point5km, out_rainfall)
message("  Saved: ", out_rainfall)


# ---- 3. LIVESTOCK DENSITY ---------------------------------------------------

message("Processing livestock density...")

livestock_rast <- rast(livestock_rast_path)
livestock_masked <- crop(livestock_rast, scotland) |> mask(scotland)
livestock_bng <- project(livestock_masked, "EPSG:27700")

template_rast <- rast(
  xmin = xmin(livestock_bng), xmax = xmax(livestock_bng),
  ymin = ymin(livestock_bng), ymax = ymax(livestock_bng),
  resolution = 2500, crs = crs(livestock_bng)
)

livestock_2point5km_rast <- resample(livestock_bng, template_rast, method = "bilinear") |>
  project("EPSG:4326")

grid_2point5km$livestock_density <- terra::extract(
  livestock_2point5km_rast, vect(grid_2point5km), fun = mean, na.rm = TRUE
)[, 2]

livestock_df <- grid_2point5km |>
  st_drop_geometry() |>
  dplyr::select(grid_id, livestock_density) |>
  mutate(
    livestock_density = as.numeric(livestock_density),
    livestock_log     = log1p(livestock_density),
    z_livestock       = scale(livestock_density)[, 1]
  )

write_csv(livestock_df, out_livestock)
message("  Saved: ", out_livestock)


# ---- 4. ELEVATION -----------------------------------------------------------

message("Processing elevation...")

elevation_rast <- rast(elevation_rast_path) |> project(crs(scotland))

grid_2point5km$elevation <- terra::extract(
  elevation_rast, vect(grid_2point5km), fun = mean, na.rm = TRUE
)[, 2]

elevation_df <- grid_2point5km |>
  st_drop_geometry() |>
  dplyr::select(grid_id, elevation) |>
  mutate(
    elevation   = as.numeric(elevation),
    z_elevation = scale(elevation)[, 1]
  )

write_csv(elevation_df, out_elevation)
message("  Saved: ", out_elevation)


# ---- 5. LANDCOVER -----------------------------------------------------------

message("Processing landcover data...")

get_landcover_area_per_class <- function(grid_sf, landcover_raster_path,
                                         num_cores = 10, chunk_size = 2000,
                                         output_dir = "temp_results") {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  base_raster <- rast(landcover_raster_path)
  
  if (!identical(crs(grid_sf), crs(base_raster))) {
    message("  Reprojecting grid to match raster CRS...")
    grid_sf <- st_transform(grid_sf, crs(base_raster))
  }
  
  # Pixel area calculation
  res_x <- res(base_raster)[1]; res_y <- res(base_raster)[2]
  if (st_crs(grid_sf)$IsGeographic) {
    center_lat <- mean(c(ext(base_raster)$ymin, ext(base_raster)$ymax))
    pixel_area_km2 <- (res_x * 111320 * cos(center_lat * pi / 180) *
                         res_y * 111320) / 1e6
  } else {
    pixel_area_km2 <- (res_x * res_y) / 1e6
  }
  message(sprintf("  Pixel area: %.6f km²", pixel_area_km2))
  
  if (num_cores > 1) plan(multisession, workers = num_cores)
  
  chunks <- split(1:nrow(grid_sf), ceiling(seq_len(nrow(grid_sf)) / chunk_size))
  
  for (i in seq_along(chunks)) {
    message(sprintf("  Chunk %d/%d", i, length(chunks)))
    chunk_results <- future_lapply(chunks[[i]], function(j) {
      tryCatch({
        raster_local <- rast(landcover_raster_path)
        extracted <- exact_extract(raster_local, grid_sf[j, ], progress = FALSE)[[1]]
        if (is.null(extracted) || nrow(extracted) == 0)
          return(tibble(grid_id = grid_sf$grid_id[j], class = NA_integer_, area_km2 = 0))
        
        val_col <- setdiff(names(extracted), "coverage_fraction")[1]
        extracted |>
          filter(!is.na(.data[[val_col]])) |>
          group_by(class = .data[[val_col]]) |>
          summarise(area_km2 = sum(coverage_fraction, na.rm = TRUE) * pixel_area_km2,
                    .groups = "drop") |>
          mutate(grid_id = grid_sf$grid_id[j])
      }, error = function(e) {
        message(sprintf("  Error in cell %d: %s", j, e$message)); NULL
      })
    }, future.seed = TRUE)
    saveRDS(bind_rows(chunk_results), file.path(output_dir, sprintf("chunk_%03d.rds", i)))
    rm(chunk_results); gc()
  }
  
  result_files <- list.files(output_dir, pattern = "^chunk_.*\\.rds$", full.names = TRUE)
  all_results  <- bind_rows(lapply(result_files, readRDS))
  unlink(output_dir, recursive = TRUE)
  
  all_results |>
    mutate(class_name = paste0("class_", class, "_km2")) |>
    dplyr::select(grid_id, class_name, area_km2) |>
    pivot_wider(names_from = class_name, values_from = area_km2, values_fill = 0)
}

landcover_areas <- get_landcover_area_per_class(
  grid_sf              = grid_2point5km,
  landcover_raster_path = landcover_rast_path,
  num_cores            = landcover_num_cores,
  chunk_size           = landcover_chunk_size
)

# Rename columns using lookup
names(landcover_areas) <- vapply(names(landcover_areas), function(col) {
  if (grepl("^class_\\d+_km2$", col)) {
    cls <- sub("^class_(\\d+)_km2$", "\\1", col)
    lbl <- landcover_lookup[[cls]]
    if (!is.na(lbl)) return(paste0(lbl, "_km2"))
  }
  col
}, character(1))

write_csv(landcover_areas, out_landcover)

# Group into ecological categories
landcover_grouped <- landcover_areas |>
  mutate(
    wetland_km2           = Bog_km2 + Fen_km2 + Saltmarsh_km2,
    woodland_km2          = Deciduous_woodland_km2 + Coniferous_woodland_km2,
    saltwater_km2         = Saltwater_km2,
    freshwater_km2        = Freshwater_km2,
    grassland_heather_km2 = Heather_grassland_km2 + Heather_km2 +
      Improved_grassland_km2 + Acid_grassland_km2 +
      Neutral_grassland_km2 + Calcareous_grassland_km2,
    urban_km2             = Urban_km2,
    suburban_km2          = Suburban_km2,
    arable_km2            = Arable_km2
  ) |>
  dplyr::select(grid_id, wetland_km2, woodland_km2, saltwater_km2, freshwater_km2,
                grassland_heather_km2, urban_km2, suburban_km2, arable_km2)

write_csv(landcover_grouped, out_landcover_grp)
message("  Saved: ", out_landcover_grp)


# ---- 6. BUILDINGS -----------------------------------------------------------

message("Processing building density...")

buildings <- vect(buildings_path) |> st_as_sf()
grid_buildings <- st_transform(grid_2point5km, crs = st_crs(buildings))

buildings_counts <- st_join(buildings, grid_buildings, join = st_within) |>
  st_drop_geometry() |>
  count(grid_id, name = "buildings_count")

buildings_df <- grid_buildings |>
  st_drop_geometry() |>
  left_join(buildings_counts, by = "grid_id") |>
  mutate(
    buildings_count = replace_na(buildings_count, 0),
    buildings_log   = log1p(buildings_count)
  ) |>
  dplyr::select(grid_id, buildings_count, buildings_log)

write_csv(buildings_df, out_buildings)
message("  Saved: ", out_buildings)


# ---- 7. POINTS OF INTEREST --------------------------------------------------

message("Processing points of interest...")

poi <- st_read(poi_path)
grid_poi <- st_transform(grid_2point5km, crs = st_crs(poi))

poi_counts <- st_join(poi, grid_poi, join = st_within) |>
  st_drop_geometry() |>
  count(grid_id, name = "poi_count")

poi_df <- grid_poi |>
  st_drop_geometry() |>
  left_join(poi_counts, by = "grid_id") |>
  mutate(
    poi_count = replace_na(poi_count, 0),
    poi_log   = log1p(poi_count)
  ) |>
  dplyr::select(grid_id, poi_count, poi_log)

write_csv(poi_df, out_poi)
message("  Saved: ", out_poi)


# ---- 8. CITIZEN SCIENCE REPORTING DENSITY -----------------------------------

message("Processing citizen science reporting density...")

cs_all <- read_csv(cs_path) |>
  mutate(
    Date_found = as.Date(Date_found, format = "%d/%m/%Y"),
    presence   = 1L
  ) |>
  drop_na(Longitude, Latitude)

cs_all_sf <- st_as_sf(cs_all, coords = c("Longitude", "Latitude"), crs = 4326) |>
  st_transform(crs = st_crs(grid_2point5km))

grid_reports <- st_join(grid_2point5km, cs_all_sf, join = st_contains) |>
  st_drop_geometry() |>
  group_by(grid_id) |>
  summarise(reports = sum(presence, na.rm = TRUE), .groups = "drop")

reports_df <- grid_2point5km |>
  st_drop_geometry() |>
  dplyr::select(grid_id) |>
  left_join(grid_reports, by = "grid_id") |>
  mutate(reports = replace_na(reports, 0))

message(sprintf("  Total reports: %d across %d grid cells",
                sum(reports_df$reports), sum(reports_df$reports > 0)))

write_csv(reports_df, out_reports)
message("  Saved: ", out_reports)


# ---- 9. BUILD FINAL COVARIATE DATASETS --------------------------------------

### SO AS TO NOT RUN THE ENTIRE SCRIPT AGAIN JUST TO REBUILD THE FINAL 
### COVARIATE DATASETS
load("01_data/covariates/dist_matrix_2point5km_meters.Rdata")
temp_2point5km <- read_csv("01_data/covariates/mean_temp_21day_2point5km.csv")
rainfall_2point5km <- read_csv("01_data/covariates/rainfall_28day_2point5km_new.csv")
livestock_df <- read_csv( "01_data/covariates/livestockDensity_2point5km_new.csv")
elevation_df <- read_csv( "01_data/covariates/elevation_2point5km_new.csv")
# landcover_areas <- read_csv("01_data/covariates/landcover_2point5km_new.csv")
landcover_grouped <- read_csv("01_data/covariates/landcover_grouped_2point5km.csv")
buildings_df <- read_csv("01_data/covariates/buildings_2point5km.csv")
poi_df <- read_csv("01_data/covariates/poi_grouped_2point5km.csv")
reports_df <- read_csv("01_data/covariates/cs_reports_2point5km.csv")

message("Building final covariate datasets...")

# --- Static land covariates
static_df <- landcover_grouped |>
  left_join(livestock_df  |> dplyr::select(grid_id, livestock_density, livestock_log), by = "grid_id") |>
  left_join(elevation_df  |> dplyr::select(grid_id, elevation),                        by = "grid_id") |>
  left_join(buildings_df  |> dplyr::select(grid_id, buildings_count, buildings_log),   by = "grid_id") |>
  left_join(poi_df        |> dplyr::select(grid_id, poi_count, poi_log),               by = "grid_id") |>
  left_join(reports_df    |> dplyr::select(grid_id, reports),                          by = "grid_id") |>
  arrange(grid_id)

z_land_df <- static_df |>
  mutate(across(
    where(is.numeric) & !matches("grid_id|area"),
    ~ scale(.x)[, 1],
    .names = "z_{.col}"
  ))

write_csv(z_land_df, out_z_land)
message("  Saved: ", out_z_land)

# --- Time-varying climate covariates
climate_df <- temp_2point5km |>
  dplyr::select(grid_id, date, min_temp_21d_celsius) |>
  left_join(rainfall_2point5km |> dplyr::select(grid_id, date, rainfall_28d),
            by = c("grid_id", "date")) |>
  group_by(grid_id) |>
  mutate(
    z_temp_min = scale(min_temp_21d_celsius)[, 1],
    z_rain     = scale(rainfall_28d)[, 1]
  ) |>
  ungroup() |>
  dplyr::select(grid_id, date, z_temp_min, z_rain)

write_csv(climate_df, out_z_climate)
message("  Saved: ", out_z_climate)

message("All covariates processed successfully!")