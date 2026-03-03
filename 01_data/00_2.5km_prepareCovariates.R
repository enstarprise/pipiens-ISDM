# ============================================================================
# SPATIAL COVARIATE EXTRACTION FOR ECOLOGICAL MODELING
# Processes climate, landcover, and human activity data at 2.5km^2 grid 
# resolution;
# tempearture and rainfall previously extracted in 01_data/gridExtractions/
# ============================================================================

# --- SETUP -------------------------------------------------------------------
library(sf)
library(terra)
library(readr)
library(dplyr)
library(exactextractr)
library(tidyr)
library(future.apply)
library(ggplot2)


#####  read csvs and distance matrix 
load("01_data/covariates/dist_matrix_2point5km_meters.Rdata")
temp_2point5km <- read_csv("01_data/covariates/min_temp_21day_2point5km_new.csv")
rainfall_2point5km <- read_csv("01_data/covariates/rainfall_28day_2point5km_new.csv")
livestock_df <- read_csv( "01_data/covariates/livestockDensity_2point5km_new.csv")
elevation_df <- read_csv( "01_data/covariates/elevation_2point5km_new.csv")
# landcover_areas <- read_csv("01_data/covariates/landcover_2point5km_new.csv")
landcover_grouped <- read_csv("01_data/covariates/landcover_grouped_2point5km.csv")
buildings_df <- read_csv("01_data/covariates/buildings_2point5km.csv")
poi_df <- read_csv("01_data/covariates/poi_grouped_2point5km.csv")
reports_df <- read_csv("01_data/covariates/cs_reports_2point5km.csv")


# --- LOAD BASE GRIDS ---------------------------------------------------------
grid_1km <- st_read("01_data/grids/grid_clipped_1km.gpkg")
grid_2point5km <- st_read("01_data/grids/grid_clipped_2point5km.gpkg")
scotland <- st_read("01_data/gridExtractions/scotlandBoundary/scotland_boundary.shp") %>% vect()

# Compute area correctly (using BNG)
grid_2point5km$area <- as.numeric(st_area(grid_2point5km)) / 1e6  # km²

# Ensure consistent CRS
# grid_2point5km <- st_transform(grid_2point5km, st_crs(scotland))
# grid_1km <- st_transform(grid_1km, st_crs(scotland))

### CREATE DIST MATRIX =======================================================
# grid centroid
grid_2point5km$centroid <- st_centroid(grid_2point5km$geom)

# Convert to projected CRS if not already (BNG should be fine)
# Get centroids as sf points
centroids_2point5km <- st_centroid(grid_2point5km)
coords_2point5km <- st_coordinates(centroids_2point5km)
# centroids_1km <- st_centroid(grid_1km)

# Create distance matrix in meters
dist_matrix_2point5km_meters <- st_distance(centroids_2point5km)
dist_matrix_2point5km_meters <- as.matrix(dist_matrix_2point5km_meters)
save(dist_matrix_2point5km_meters, file = "dist_matrix_2point5km_meters.Rdata")

summary(as.vector(dist_matrix_2point5km_meters)) 
max(dist_matrix_2point5km_meters) 

# Convert to matrix and km if needed
# dist_matrix_2point5km_km <- as.matrix(dist_matrix_2point5km_meters) / 1000  # to km
# save(dist_matrix_2point5km_km, file = "dist_matrix_2point5km_km.Rdata")


## PREPARE GRID 2KM (CLIPPED) AS A WKT CSV FILE FOR USE IN THE CLUSTER
grid_2point5km_wkt <- grid_2point5km |>
  mutate(geometry = st_as_text(geom))  # convert geometry to readable WKT

write_csv(grid_2point5km_wkt, "01_data/grids/grid_clipped_2point5km_wkt.csv")
write_csv(grid_2point5km_wkt, "01_data/csvs/grid_clipped_2point5km_wkt.csv")

# --- 1. CLIMATE DATA: TEMPERATURE --------------------------------------------
message("Processing minimum temperature data...")

temp_1km <- read_csv("01_data/csvs/tasmin_grid_21day.csv")

# 1. Map 1km grids to 2.5km grids
mapping_1km_to_2point5km <- st_join(
  grid_1km %>% dplyr::select(grid_id), 
  grid_2point5km %>% dplyr::select(grid_id), 
  left = FALSE
) %>%
  st_drop_geometry() %>%
  rename(
    grid_id_1km = grid_id.x,        # 1km grid
    grid_id_2point5km = grid_id.y   # 2.5km grid
  )

# 2. Aggregate temperature to 2.5km grids
temp_2point5km <- temp_1km %>%
  # Join the mapping to get the 2.5km grid ID
  left_join(mapping_1km_to_2point5km, by = c("grid_id" = "grid_id_1km")) %>%
  
  # Group by 2.5km grid and date
  group_by(grid_id_2point5km, date) %>%
  
  # Compute mean of 7-day min temperature
  summarise(
    min_temp_21d_celsius = mean(mean_min_temp_7d_celsius, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  
  #standardize the temperature
  mutate(
    min_temp_21d_celsius = as.numeric(min_temp_21d_celsius),
    z_temp_min = scale(min_temp_21d_celsius)[,1]
  )

names(temp_2point5km)[1] <- "grid_id"


write_csv(temp_2point5km, "01_data/covariates/min_temp_21day_2point5km_new.csv")


# --- 2. CLIMATE DATA: RAINFALL -----------------------------------------------
message("Processing rainfall data...")

rainfall_1km <- read_csv("01_data/csvs/rainfall_grid_new.csv")

# Map and aggregate
rainfall_2km <- rainfall_1km %>%
  left_join(mapping_1km_to_2km, by = "grid_id") %>%
  group_by(grid_id_area2km, date) %>%
  summarise(rainfall_28d = mean(rainfall_28d, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    rainfall_28d = as.numeric(rainfall_28d),
    z_rain = scale(rainfall_28d)[,1]
  )

rainfall_2point5km <- rainfall_1km %>%
  # Join the mapping to get the 2.5km grid ID
  left_join(mapping_1km_to_2point5km, by = c("grid_id" = "grid_id_1km")) %>%
  
  # Group by 2.5km grid and date
  group_by(grid_id_2point5km, date) %>%
  
  # Compute mean of 21-day mean rainfall
  summarise(rainfall_28d = mean(rainfall_28d, na.rm = TRUE), .groups = "drop") %>%
  
  #standardize the temperature
  mutate(
    rainfall_28d = as.numeric(rainfall_28d),
    z_rain = scale(rainfall_28d)[,1]
  )

names(rainfall_2point5km)[1] <- "grid_id"

write_csv(rainfall_2point5km, "01_data/covariates/rainfall_28day_2point5km_new.csv")


# --- 3. LIVESTOCK DENSITY ----------------------------------------------------
message("Processing livestock density...")

livestock_rast <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/livestockDensity/GLW4-2020.D-DA.GLEAM3-ALL-LU.tif")

# 1. Crop and mask livestock raster (still in WGS84)
livestock_masked <- crop(livestock_rast, scotland) %>% mask(scotland)

# 2. Reproject to British National Grid (meters)
livestock_bng <- project(livestock_masked, "EPSG:27700")

# 3. Create a template raster at 2.5 km resolution
res_target <- 2500  # meters
template_rast <- rast(
  xmin = xmin(livestock_bng),
  xmax = xmax(livestock_bng),
  ymin = ymin(livestock_bng),
  ymax = ymax(livestock_bng),
  resolution = res_target,
  crs = crs(livestock_bng)
)

# 4. Resample livestock raster to 2.5 km cells
livestock_2point5km_bng <- resample(livestock_bng, template_rast, method = "bilinear")

# 5. (Optional) Reproject back to WGS84
livestock_2point5km <- project(livestock_2point5km_bng, "EPSG:4326")

# 6. Extract livestock values to your grid
grid_2point5km$livestock_density <- terra::extract(
  livestock_2point5km, 
  vect(grid_2point5km), 
  fun = mean, 
  na.rm = TRUE
)[, 2]

# 7. Convert to a clean dataframe
livestock_df <- grid_2point5km %>%
  st_drop_geometry() %>%
  dplyr::select(grid_id, livestock_density) %>%
  mutate(
    livestock_density = as.numeric(livestock_density),
    livestock_log = log1p(livestock_density),   # log(1 + x)
    z_livestock = scale(livestock_density)[,1]
  )

write_csv(livestock_df, "01_data/covariates/livestockDensity_2point5km_new.csv")


# --- 4. ELEVATION ------------------------------------------------------------
message("Processing elevation...")

elevation <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/elevation_DEM/scotland_elevation.tif")
elevation <- project(elevation, crs(scotland)) # 4362

grid_2point5km$elevation <- terra::extract(
  elevation, 
  vect(grid_2point5km), 
  fun = mean, 
  na.rm = TRUE
)[, 2]

elevation_df <- grid_2point5km %>%
  st_drop_geometry() %>%
  dplyr::select(grid_id, elevation) %>%
  mutate(
    elevation = as.numeric(elevation),
    z_elevation = scale(elevation)[,1]
  )

write_csv(elevation_df, "01_data/covariates/elevation_area2km_2point5km_new.csv")


# --- 5. LANDCOVER CLASSIFICATION ---------------------------------------------
message("Processing landcover data...")

# Helper function for accurate area calculation
get_landcover_area_per_class <- function(
    grid_sf,
    landcover_raster_path,
    num_cores = 10,
    chunk_size = 2000,
    output_dir = "temp_results"
) {
  
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  base_raster <- rast(landcover_raster_path)
  
  if (!identical(crs(grid_sf), crs(base_raster))) {
    message("Reprojecting grid to match raster CRS...")
    grid_sf <- st_transform(grid_sf, crs(base_raster))
  }
  
  # Calculate pixel area
  res_x <- res(base_raster)[1]
  res_y <- res(base_raster)[2]
  crs_info <- st_crs(grid_sf)
  
  if (crs_info$IsGeographic) {
    ext_bbox <- ext(base_raster)
    center_lat <- (ext_bbox$ymin + ext_bbox$ymax) / 2
    lat_m_per_deg <- 111320
    lon_m_per_deg <- 111320 * cos(center_lat * pi / 180)
    pixel_area_km2 <- (res_x * lon_m_per_deg * res_y * lat_m_per_deg) / 1e6
  } else {
    pixel_area_km2 <- (res_x * res_y) / 1e6
  }
  
  message(sprintf("Calculated pixel area: %.6f km²", pixel_area_km2))
  
  # Parallel processing setup
  if (num_cores > 1) {
    plan(multisession, workers = num_cores)
  }
  
  chunks <- split(1:nrow(grid_sf), ceiling(seq_along(1:nrow(grid_sf)) / chunk_size))
  
  # Process chunks
  for (i in seq_along(chunks)) {
    message(sprintf("Processing chunk %d/%d", i, length(chunks)))
    
    chunk_results <- future_lapply(
      chunks[[i]],
      function(j) {
        tryCatch({
          raster_local <- rast(landcover_raster_path)
          extracted <- exact_extract(raster_local, grid_sf[j, ], progress = FALSE)[[1]]
          
          if (is.null(extracted) || nrow(extracted) == 0) {
          
            return(tibble(grid_id = grid_sf$grid_id[j], class = NA_integer_, area_km2 = 0))
          }
          
          value_col <- setdiff(names(extracted), "coverage_fraction")[1]
          
          extracted %>%
            filter(!is.na(.data[[value_col]])) %>%
            group_by(class = .data[[value_col]]) %>%
            summarise(area_km2 = sum(coverage_fraction, na.rm = TRUE) * pixel_area_km2, .groups = "drop") %>%
            mutate(grid_id = grid_sf$grid_id[j])  # FIX: Use actual ID
        }, error = function(e) {
          message(sprintf("Error in cell %d: %s", j, e$message))
          NULL
        })
      },
      future.seed = TRUE
    )
    
    saveRDS(bind_rows(chunk_results), file.path(output_dir, sprintf("chunk_%03d.rds", i)))
    rm(chunk_results)
    gc()
  }
  
  # Combine results
  result_files <- list.files(output_dir, pattern = "^chunk_.*\\.rds$", full.names = TRUE)
  all_results <- bind_rows(lapply(result_files, readRDS))
  unlink(output_dir, recursive = TRUE)
  
  # Pivot to wide format
  all_results %>%
    mutate(class_name = paste0("class_", class, "_km2")) %>%
    dplyr::select(grid_id, class_name, area_km2) %>%  # FIX: Use grid_id
    pivot_wider(names_from = class_name, values_from = area_km2, values_fill = 0)
}

# Process landcover
landcover <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/landCover/UK_CEH/data/7727ce7d-531e-4d77-b756-5cc59ff016bd/ukregion-scotland.tif")
# writeRaster(landcover, "01_data/gridExtractions/landcover_27700.tif", overwrite = TRUE)

landcover_areas <- get_landcover_area_per_class(
  grid_sf = grid_2point5km,  
  landcover_raster_path = "01_data/gridExtractions/landcover_27700.tif",
  num_cores = 10,
  chunk_size = 2000
)

# Rename columns using lookup
landcover_lookup <- c(
  `0`  = "Unnamed",`1`  = "Deciduous_woodland",`2`  = "Coniferous_woodland",
  `3`  = "Arable",`4`  = "Improved_grassland",`5`  = "Neutral_grassland",
  `6`  = "Calcareous_grassland",`7`  = "Acid_grassland",`8`  = "Fen",
  `9`  = "Heather",`10` = "Heather_grassland",`11` = "Bog",`12` = "Inland_rock",
  `13` = "Saltwater",`14` = "Freshwater",`15` = "Supralittoral_rock",
  `16` = "Supralittoral_sediment",`17` = "Littoral_rock",`18` = "Littoral_sediment",
  `19` = "Saltmarsh",`20` = "Urban",`21` = "Suburban"
)

# Rename columns using the lookup
names(landcover_areas) <- vapply(
  names(landcover_areas),
  function(col) {
    if (grepl("^class_\\d+_km2$", col)) {
      class_num <- sub("^class_(\\d+)_km2$", "\\1", col)
      class_label <- landcover_lookup[[class_num]]
      if (!is.na(class_label)) {
        return(paste0(class_label, "_km2"))
      }
    }
    return(col)
  },
  character(1)
)

write_csv(landcover_areas, "01_data/covariates/landcover_2point5km_new.csv")

# Group landcover into ecological categories
landcover_grouped <- landcover_areas %>%
  mutate(
    wetland_km2 = Bog_km2 + Fen_km2 + Saltmarsh_km2,
    woodland_km2 = Deciduous_woodland_km2 + Coniferous_woodland_km2,
    saltwater_km2 = Saltwater_km2,
    freshwater_km2 = Freshwater_km2,
    grassland_heather_km2 = Heather_grassland_km2 + Heather_km2 + 
      Improved_grassland_km2 + Acid_grassland_km2 + 
      Neutral_grassland_km2 + Calcareous_grassland_km2,
    urban_km2 = Urban_km2, 
    suburban_km2 = Suburban_km2,
    arable_km2 = Arable_km2
  ) %>%
  rename(grid_id = grid_id) %>%  # Standardize column name
  dplyr::select(grid_id, wetland_km2, woodland_km2, saltwater_km2, freshwater_km2, 
                grassland_heather_km2, urban_km2, suburban_km2, arable_km2)


### log transform wetland
# landcover_grouped <- landcover_grouped %>%
#   mutate(freshwater_km2_log = log1p(freshwater_km2))   # log(1 + x))


write_csv(landcover_grouped, "01_data/covariates/landcover_grouped_2point5km.csv")


# --- 8. BUILDING DENSITY -----------------------------------------------------
message("Processing building density data...")
buildings <- vect("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/scotlandOSMBuildings_Others.shp/gis_osm_buildings_a_free_1.shp") %>%
  st_as_sf()


# Count buildings per grid cell (using full poi dataset)
grid_2point5km_buildings <- st_transform(grid_2point5km, crs = st_crs(buildings))
buildings_counts <- st_join(buildings, grid_2point5km_buildings, join = st_within) %>%
  st_drop_geometry() %>%
  count(grid_id, name = "buildings_count")

buildings_df <- grid_2point5km_buildings %>%
  st_drop_geometry() %>%
  left_join(buildings_counts, by = "grid_id") %>%
  mutate(
    buildings_count = replace_na(buildings_count, 0),
    buildings_log = log1p(buildings_count),   # log(1 + x)
    # buildings_count_grouped = pmin(buildings_count, 100)  # Cap at 100
  ) %>%
  dplyr::select(grid_id, buildings_count, buildings_log)

write_csv(buildings_df, "01_data/covariates/buildings_2point5km.csv")


# --- 7. POINTS OF INTEREST (POI) ---------------------------------------------
message("Processing points of interest...")

poi <- st_read("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/POI/poi_uk.gpkg")

# Filter relevant categories
# relevant_categories <- c(
#   "park", "beach", "trail", "nature_reserve", "hiking_trail", 
#   "campground", "lake", "river", "forest", "botanical_garden",
#   "national_park", "state_park", "playground", "dog_park", 
#   "golf_course", "restaurant", "cafe", "pub", "hotel"
# )
# 
# poi_subset <- poi %>% filter(main_category %in% relevant_categories)

# Count POIs per grid cell (using full poi dataset)
grid_2point5km_poi <- st_transform(grid_2point5km, crs = st_crs(poi))
poi_counts <- st_join(poi, grid_2point5km_poi, join = st_within) %>%
  st_drop_geometry() %>%
  count(grid_id, name = "poi_count")

poi_df <- grid_2point5km_poi %>%
  st_drop_geometry() %>%
  left_join(poi_counts, by = "grid_id") %>%
  mutate(
    poi_count = replace_na(poi_count, 0),
    poi_log = log1p(poi_count),   # log(1 + x)
    # poi_count_grouped = pmin(poi_count, 100)  # Cap at 100
  ) %>%
  dplyr::select(grid_id, poi_count, poi_log)

write_csv(poi_df, "01_data/covariates/poi_grouped_2point5km.csv")

# --- 8. TOTAL REPORTING ACTIVITY?INTENSITY ------------------------------------
message("Processing citizen science reporting density...")

# Load all citizen science reports (not just presences)
cs_all <- read_csv("01_data/csvs/cs_df.csv")

# Clean dates and add temporal variables
cs_all <- cs_all %>%
  mutate(
    Date_found = as.Date(Date_found, format = "%d/%m/%Y"),
    Year = format(Date_found, "%Y"),
    Month = format(Date_found, "%m"),
    presence = 1  # Each observation = 1 report
  ) %>%
  drop_na(Longitude, Latitude)  # Remove records with missing coordinates

# Convert to spatial object
cs_all_sf <- st_as_sf(cs_all, coords = c("Longitude", "Latitude"), crs = 4326)
cs_all_sf <- st_transform(cs_all_sf, crs = st_crs(grid_2km))

# Count total reports per grid cell
grid_reports <- st_join(grid_2point5km, cs_all_sf, join = st_contains) %>%
  st_drop_geometry() %>%
  group_by(grid_id) %>%
  summarise(reports = sum(presence, na.rm = TRUE), .groups = "drop")

# Join back to grid and fill missing with 0
reports_df <- grid_2point5km %>%
  st_drop_geometry() %>%
  dplyr::select(grid_id) %>%
  left_join(grid_reports, by = "grid_id") %>%
  mutate(reports = replace_na(reports, 0))

message(sprintf("  Total reports: %d across %d grid cells", 
                sum(reports_df$reports), 
                sum(reports_df$reports > 0)))

write_csv(reports_df, "01_data/covariates/cs_reports_2point5km.csv")

# --- 9. BUILD FINAL COVARIATE DATASETS --------------------------------------
message("Building final covariate datasets...")

build_covariates <- function(
    dist_matrix,
    landcover_grouped,
    livestock_df,
    elevation_df,
    buildings_df,
    poi_df,
    reports_df,
    temp_df,
    rain_df,
    land_out = "01_data/covariates/z_land_2point5km.csv",
    climate_out = "01_data/covariates/z_climate_2point5km.csv"
) {
  
  # --- Distance matrix ---
  message("Distance matrix...")
  # save(dist_matrix_2point5km_meters, file = "dist_matrix_2point5km_meters.Rdata")
  
  
  # --- Static land covariates ---
  message("Building static land covariates...")
  
  static_df <- landcover_grouped %>%
    left_join(livestock_df %>% dplyr::select(grid_id, livestock_density), by = "grid_id") %>%
    left_join(livestock_df %>% dplyr::select(grid_id, livestock_log), by = "grid_id") %>% ## z score the log column
    left_join(elevation_df %>% dplyr::select(grid_id, elevation), by = "grid_id") %>%
    left_join(buildings_df %>% dplyr::select(grid_id, buildings_count), by = "grid_id") %>%
    left_join(poi_df %>% dplyr::select(grid_id, poi_count), by = "grid_id") %>%
    left_join(buildings_df %>% dplyr::select(grid_id, buildings_log), by = "grid_id") %>% ## z score the log column
    left_join(poi_df %>% dplyr::select(grid_id, poi_log), by = "grid_id") %>% ## z score the log column
    left_join(reports_df %>% dplyr::select(grid_id, reports), by = "grid_id") %>%
    arrange(grid_id)
  
  # Auto z-score all numeric columns EXCEPT grid_id and area columns
  z_land_df <- static_df %>%
    mutate(
      across(
        where(is.numeric) & !matches("grid_id|area"),  # Exclude grid_id and area from scaling
        ~scale(.x)[,1], 
        .names = "z_{.col}"
      )
    )
  
  # also z-score the wetland column?
  
  write_csv(z_land_df, land_out)
  message(paste("Wrote:", land_out))
  
  # --- Time-varying climate covariates ---
  message("Building time-varying climate covariates...")
  
  climate_df <- temp_df %>%
    dplyr::select(grid_id, date, min_temp_21d_celsius) %>%
    left_join(rain_df %>% dplyr::select(grid_id, date, rainfall_28d), 
              by = c("grid_id", "date")) %>%
    group_by(grid_id) %>%
    mutate(
      z_temp_min = scale(min_temp_21d_celsius)[,1],
      z_rain = scale(rainfall_28d)[,1]
    ) %>%
    ungroup() %>%
    dplyr::select(grid_id, date, z_temp_min, z_rain)
  
  write_csv(climate_df, climate_out)
  message(paste("Wrote:", climate_out))
  
  list(z_land = z_land_df, z_climate = climate_df)
}

# Execute final build
final_covariates <- build_covariates(
  dist_matrix = dist_matrix_2point5km_meters,
  landcover_grouped = landcover_grouped,
  livestock_df = livestock_df,
  elevation_df = elevation_df,
  buildings_df = buildings_df,
  poi_df = poi_df,
  reports_df = reports_df,
  temp_df = temp_2point5km,
  rain_df = rainfall_2point5km
)

message("✓ All covariates processed successfully!")



