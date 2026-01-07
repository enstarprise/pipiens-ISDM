library(sf)
library(raster)
library(terra)
library(readr)
library(exactextractr)
library(dplyr)
grid_clipped <- st_read("grids/grid_clipped_area2km.gpkg")

scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp") %>%
  vect()

grid_clipped <- st_transform(grid_clipped, st_crs(scotland))
st_crs(grid_clipped)


grid_clipped <- vect(grid_clipped)
grid_clipped <- st_as_sf(grid_clipped)

grid_1km <- st_read("grids/grid_clipped_1km.gpkg")
st_crs(grid_clipped) == st_crs(grid_1km)
grid_1km <- st_transform(grid_1km, st_crs(scotland))

### --- mean temperature---
library(dplyr)
library(sf)

# grid_1km_temp <- read_csv("csvs/tasmean_grid_new.csv")
# 
# grid_2point5 <- grid_clipped %>%
#   rename(grid_id_2point5km = grid_id)
# 
# mapping_1km_to_2_5km <- st_join(
#   grid_1km %>% select(grid_id), 
#   grid_2point5 %>% select(grid_id_2point5km), 
#   left = FALSE
# ) %>%
#   st_drop_geometry()
# 
# temp_df_mapped <- grid_1km_temp %>%
#   left_join(mapping_1km_to_2_5km, by = "grid_id")
# 
# temp_2_5km <- temp_df_mapped %>%
#   group_by(grid_id_2point5km, date) %>%
#   summarise(mean_temp_7d_celsius = mean(mean_temp_7d_celsius, na.rm = TRUE),
#             .groups = "drop")
# 
# length(unique(temp_2_5km$grid_id_2point5km))

#### MINIMUM TEMPERATURE IINSTEDAD OF MEAN
grid_1km_temp <- read_csv("01_data/csvs/tasmin_grid_21day.csv")

grid_area2km <- grid_clipped %>%
  rename(grid_id_area2km = grid_id)

mapping_1km_to_2km <- st_join(
  grid_1km %>% dplyr::select(grid_id), 
  grid_area2km %>% dplyr::select(grid_id_area2km), 
  left = FALSE
) %>%
  st_drop_geometry()

temp_df_mapped <- grid_1km_temp %>%
  left_join(mapping_1km_to_2km, by = "grid_id")

temp_2km <- temp_df_mapped %>%
  group_by(grid_id_area2km, date) %>%
  summarise(min_temp_21d_celsius = mean(mean_min_temp_7d_celsius, na.rm = TRUE),
            .groups = "drop")

length(unique(temp_2km$grid_id_area2km))

temp_2km <- temp_2km %>%
  mutate(min_temp_21d_celsius = as.numeric(min_temp_21d_celsius)) %>%
  mutate(z_temp_min = scale(min_temp_21d_celsius)[,1])


write.csv(temp_2km, "01_data/covariates/min_temp_21day_area2km.csv")

### --- cumulative rainfall ---
rainfall_1km <- read_csv("01_data/csvs/rainfall_grid_new.csv")

# grid_1km <- st_transform(grid_1km, st_crs(grid_2point5))

rainfall_mapping_1km_to_2km <- st_join(
  grid_1km %>% dplyr::select(grid_id),
  grid_area2km %>% dplyr::select(grid_id_area2km),
  join = st_intersects
) %>%
  st_drop_geometry()

rainfall_mapped <- rainfall_1km %>%
  left_join(rainfall_mapping_1km_to_2km, by = "grid_id")
length(unique(rainfall_mapped$grid_id_2km))


rainfall_area2km <- rainfall_mapped %>%
  group_by(grid_id_area2km, date) %>%
  summarise(
    rainfall_28d = mean(rainfall_28d, na.rm = TRUE),
    .groups = "drop"
  )

length(unique(rainfall_area2km$grid_id_area2km))

rainfall_area2km <- rainfall_area2km %>%
  mutate(rainfall_28d = as.numeric(rainfall_28d)) %>%
  mutate(z_rain = scale(rainfall_28d)[,1])

write.csv(rainfall_area2km, "01_data/covariates/rainfall_28day_area2km.csv")

# merge rainfall and temperature
# climate_data_2point5km <- full_join(rainfall_2_5km, temp_2_5km, by = c("grid_id_2point5km", "date"))
# length(unique(climate_data_2point5km$grid_id_2point5km))
# 
# write.csv(climate_data_2point5km ,"csvs/climate_data_2point5km.csv")

### --- human popualtion density ---
# humanDensity <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R3_C18/GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0_R3_C18.tif")
# humanDensity <- project(humanDensity, scotland)
# st_crs(humanDensity) == st_crs(grid_clipped)
# human_density_2_5km <- exact_extract(humanDensity, 
#                                      grid_clipped, fun = "mean", 
#                                      progress = TRUE)
# 
# grid_clipped$human_density <- human_density_2_5km


### --- livestock density ---
livestockDensity <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/livestockDensity/GLW4-2020.D-DA.GLEAM3-ALL-LU.tif")
crs(livestockDensity) <- "EPSG:4326" # or whatever its original CRS is
st_crs(livestockDensity) == st_crs(grid_clipped)

# crop and mask to scotland's shpaefile
livestockDensity <- crop(livestockDensity,  ext(scotland))
livestock_masked <- mask(livestockDensity, scotland)
livestockDensity <- raster(livestock_masked)

# terra format for extraction
livestock_terra <- rast(livestockDensity)  # convert to spatraster
# resample 
# 10 km / 1.4 km ≈ 7.142857
livestock_area2km <- disagg(rast(livestockDensity), fact = 7.142857, method = "bilinear")

grid_area2km_livestock <- grid_clipped

grid_area2km_livestock$livestock_density <- terra::extract(
  livestock_area2km, #livestock_1km
  grid_clipped, 
  fun = mean,
  na.rm = TRUE
)[, 2]


grid_clipped <- st_join(grid_clipped, 
                        grid_area2km_livestock[,c("grid_id", "livestock_density")])


grid_area2km_livestock <- read_csv("01_data/covariates/livestockDensity_area2km.csv")
grid_area2km_livestock <- grid_area2km_livestock %>%
  mutate(
    livestock_density = as.numeric(livestock_density)
  ) %>%
  mutate(
    z_livestock   = scale(livestock_density)[,1]
  )

z_livestock <- grid_area2km_livestock %>% distinct(grid_id, z_livestock) %>% arrange(grid_id) %>% pull(z_livestock)

write_csv(grid_area2km_livestock, "01_data/covariates/livestockDensity_area2km.csv")


### --- elevation ---
elevation <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/elevation_DEM/scotland_elevation.tif")
elevation <- project(elevation, crs(scotland))
st_crs(elevation) == st_crs(grid_clipped)
# grid_terra <- vect(grid_clipped) # convrt to spatvector; terra format for extraction

grid_area2km_elevation2 <- grid_clipped
grid_area2km_elevation2$elevation <- terra::extract(
  elevation, 
  grid_clipped, 
  fun = mean,
  na.rm = TRUE
)[, 2] # will take a while 

grid_area2km_elevation2 <- grid_area2km_elevation2 %>%
  mutate(
  elevation = as.numeric(elevation)
  ) %>%
  mutate(
    z_elevation   = scale(elevation)[,1]
    )



write_csv(grid_area2km_elevation2, "01_data/covariates/elevation_area2km.csv")
# grid_area2km_elevation2 <- read_csv("01_data/covariates/elevation_area2km.csv")

z_elevation <- grid_area2km_elevation2  %>% distinct(grid_id, z_elevation) %>% arrange(grid_id) %>% pull(z_elevation)
write_csv(z_elevation, "01_data/z_elevation.csv")

### --- landcover km2 proportions ---

# landcover raster already loaded as 'landcover'
landcover <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/landCover/UK_CEH/data/7727ce7d-531e-4d77-b756-5cc59ff016bd/ukregion-scotland.tif")
grid_clipped <- st_read("grids/grid_clipped_area2km.gpkg")
# landcover <- project(landcover, grid_clipped)
# st_crs(landcover) == st_crs(grid_clipped)


# Reproject to EPSG:27700 for exact calculations
# landcover_proj <- terra::project(landcover, "EPSG:27700", method = "near")
# grid_clipped_proj <- st_transform(grid_clipped, "EPSG:27700")

get_landcover_area_per_class_fixed <- function(
    grid_sf,
    landcover_raster_path,
    num_cores = 10,
    chunk_size = 2000,
    output_dir = "temp_results"
) {
  require(terra)
  require(sf)
  require(future.apply)
  require(exactextractr)
  require(dplyr)
  require(tidyr)
  
  if (!inherits(grid_sf, "sf")) stop("grid_sf must be an sf object")
  if (!file.exists(landcover_raster_path)) stop("Raster file path does not exist")
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  if (!"grid_id" %in% names(grid_sf)) {
    grid_sf$grid_id <- seq_len(nrow(grid_sf))
  }
  
  base_raster <- terra::rast(landcover_raster_path)
  
  if (!identical(crs(grid_sf), crs(base_raster))) {
    message("Reprojecting grid to match raster CRS...")
    grid_sf <- st_transform(grid_sf, crs(base_raster))
  }
  
  # Calculate pixel area properly
  # Get resolution in map units
  res_x <- terra::res(base_raster)[1]
  res_y <- terra::res(base_raster)[2]
  
  # Check if CRS is projected or geographic
  crs_info <- st_crs(grid_sf)
  
  if (crs_info$IsGeographic) {
    # For geographic CRS, calculate actual area accounting for latitude
    message("Geographic CRS detected. Calculating accurate pixel areas...")
    
    # Get the center latitude of the raster
    ext <- terra::ext(base_raster)
    center_lat <- (ext$ymin + ext$ymax) / 2
    
    # Convert degrees to km at this latitude
    # Latitude: always ~111.32 km per degree
    # Longitude: varies with latitude: 111.32 * cos(lat in radians)
    lat_km_per_deg <- 111.32
    lon_km_per_deg <- 111.32 * cos(center_lat * pi / 180)
    
    # For very fine resolution rasters, calculate in meters first for precision
    lat_m_per_deg <- lat_km_per_deg * 1000
    lon_m_per_deg <- lon_km_per_deg * 1000
    
    # Calculate pixel dimensions in meters
    pixel_height_m <- res_y * lat_m_per_deg
    pixel_width_m <- res_x * lon_m_per_deg
    
    # Convert to km²
    pixel_area_km2 <- (pixel_width_m * pixel_height_m) / 1e6
    
    message(sprintf("At latitude %.2f°: lon=%.2f km/deg, lat=%.2f km/deg", 
                    center_lat, lon_km_per_deg, lat_km_per_deg))
    message(sprintf("Resolution: %.6f° x %.6f° (%.1f m × %.1f m)", 
                    res_x, res_y, pixel_width_m, pixel_height_m))
  } else {
    # For projected CRS, convert map units to km²
    # Check units
    units <- crs_info$units_gdal
    
    if (is.null(units) || units == "metre" || units == "m") {
      # Resolution is in meters, convert to km²
      pixel_area_km2 <- (res_x * res_y) / 1e6
      message(sprintf("Projected CRS with meters. Pixel area: %.6f km²", pixel_area_km2))
    } else if (units == "US survey foot" || units == "foot") {
      # Convert feet to meters then to km²
      pixel_area_km2 <- (res_x * 0.3048 * res_y * 0.3048) / 1e6
      message(sprintf("Projected CRS with feet. Pixel area: %.6f km²", pixel_area_km2))
    } else {
      warning(sprintf("Unknown units: %s. Assuming meters.", units))
      pixel_area_km2 <- (res_x * res_y) / 1e6
    }
  }
  
  message(sprintf("Calculated pixel area: %.6f km²", pixel_area_km2))
  
  if (num_cores > 1) {
    future::plan(future::multisession, workers = num_cores)
  } else {
    future::plan(future::sequential)
  }
  
  chunks <- split(1:nrow(grid_sf), ceiling(seq_along(1:nrow(grid_sf)) / chunk_size))
  
  for (i in seq_along(chunks)) {
    message(sprintf("Processing chunk %d/%d", i, length(chunks)))
    
    chunk_results <- future_lapply(
      chunks[[i]],
      function(j) {
        tryCatch({
          raster_local <- terra::rast(landcover_raster_path)
          extracted <- exactextractr::exact_extract(
            raster_local, 
            grid_sf[j, ], 
            progress = FALSE
          )[[1]]
          
          if (is.null(extracted) || nrow(extracted) == 0) {
            # no overlap → return zero-area placeholder
            return(
              tibble(
                grid_id = grid_sf$grid_id[j],
                class = NA_integer_,
                area_km2 = 0
              )
            )
          }
          
          
          # Dynamically detect value column
          value_col <- setdiff(names(extracted), "coverage_fraction")[1]
          
          extracted %>%
            filter(!is.na(.data[[value_col]])) %>%
            group_by(class = .data[[value_col]]) %>%
            summarise(
              area_km2 = sum(coverage_fraction, na.rm = TRUE) * pixel_area_km2, 
              .groups = "drop"
            ) %>%
            mutate(grid_id = grid_sf$grid_id[j])
        }, error = function(e) {
          message(sprintf("Error in cell %d: %s", j, e$message))
          NULL
        })
      },
      future.seed = TRUE
    )
    
    chunk_df <- bind_rows(chunk_results)
    
    saveRDS(chunk_df, file.path(output_dir, sprintf("chunk_%03d.rds", i)))
    rm(chunk_results)
    gc()
  }
  
  result_files <- list.files(output_dir, pattern = "^chunk_.*\\.rds$", full.names = TRUE)
  all_results <- bind_rows(lapply(result_files, readRDS))
  unlink(output_dir, recursive = TRUE)
  
  wide_result <- all_results %>%
    mutate(class_name = paste0("class_", class, "_km2")) %>%
    dplyr::select(grid_id, class_name, area_km2) %>%
    tidyr::pivot_wider(
      names_from = class_name,
      values_from = area_km2,
      values_fill = 0
    )
  
  return(wide_result)
}

# Usage:
terra::writeRaster(landcover, "landcover_27700.tif", overwrite = TRUE)
wide_class_areas <- get_landcover_area_per_class_fixed(
  grid_sf = grid_clipped,
  landcover_raster_path = "landcover_27700.tif",
  num_cores = 10,
  chunk_size = 2000
)

# pixel area of 0.0001 km² = 100 m² means  pixels are 10m × 10m: 
#  10m × 10m = 100 m² = 0.0001 km² 
# 1,000,000 m² ÷ 100 m² per pixel = 10,000 pixels per grid cell

# Check the results
class_cols <- grep("^class_.*_km2$", names(wide_class_areas), value = TRUE)
wide_class_areas$total_area <- rowSums(wide_class_areas[, class_cols], na.rm = TRUE)

# Summary statistics
summary(wide_class_areas$total_area)

# Should show values close to 1.0 km² for most grid cells
hist(wide_class_areas$total_area, main = "Total Area per Grid Cell", xlab = "Area (km²)")


landcover_lookup <- c(
  `0`  = "Unnamed",
  `1`  = "Deciduous_woodland",
  `2`  = "Coniferous_woodland",
  `3`  = "Arable",
  `4`  = "Improved_grassland",
  `5`  = "Neutral_grassland",
  `6`  = "Calcareous_grassland",
  `7`  = "Acid_grassland",
  `8`  = "Fen",
  `9`  = "Heather",
  `10` = "Heather_grassland",
  `11` = "Bog",
  `12` = "Inland_rock",
  `13` = "Saltwater",
  `14` = "Freshwater",
  `15` = "Supralittoral_rock",
  `16` = "Supralittoral_sediment",
  `17` = "Littoral_rock",
  `18` = "Littoral_sediment",
  `19` = "Saltmarsh",
  `20` = "Urban",
  `21` = "Suburban"
)

# Rename columns using the lookup
names(wide_class_areas) <- vapply(
  names(wide_class_areas),
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

write.csv(wide_class_areas, "01_data/covariates/landcover_area2km.csv", row.names = FALSE)
# wide_class_areas <- read_csv("csvs/landcover_2point5km.csv")

library(ggplot2)

# Sum total area across all classes per grid cell
total_area_df <- grid_clipped %>%
  left_join(
    wide_class_areas %>%
      mutate(total_km2 = rowSums(dplyr::select(., starts_with("class_")))),
    by = "grid_id"
  )


ggplot(total_area_df) +
  geom_sf(aes(fill = total_km2)) +
  scale_fill_viridis_c() +
  theme_minimal() +
  labs(title = "Total Landcover Area per Grid Cell", fill = "Area (km²)")

extracted <- exactextractr::exact_extract(landcover, grid_clipped[1, ])
sum(extracted[[1]]$coverage_fraction) * pixel_area_km2


raster_extent <- st_as_sf(as.polygons(ext(landcover), crs = crs(landcover)))
grid_clipped <- st_intersection(grid_clipped, raster_extent)

summary(total_area_df$total_km2)

# 
# wide_class_areas <- read_csv("01_data/covariates/landcover_area2km.csv")

landcover_full <- grid_clipped %>%
  st_drop_geometry() %>%
  dplyr::select(grid_id) %>%
  left_join(landcover_grouped, by = "grid_id")


landcover_grouped <- wide_class_areas %>%
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
  dplyr::select(-Deciduous_woodland_km2, -Coniferous_woodland_km2,
                -Supralittoral_rock_km2, -Supralittoral_sediment_km2,
                -Heather_grassland_km2, -Heather_km2, -Improved_grassland_km2,
                -Acid_grassland_km2, -Neutral_grassland_km2, -Calcareous_grassland_km2,
                -Littoral_rock_km2, -Littoral_sediment_km2, -Urban_km2, -Suburban_km2, -Arable_km2, 
                -Inland_rock_km2, -Bog_km2, -Fen_km2, -Saltmarsh_km2, -Saltwater_km2, -Freshwater_km2,
                -Unnamed_km2) %>%
  mutate(
    wetland_km2   = as.numeric(wetland_km2),
    saltwater_km2 = as.numeric(saltwater_km2),
    freshwater_km2 = as.numeric(freshwater_km2),
    woodland_km2 = as.numeric(woodland_km2),
    grassland_heather_km2 = as.numeric(grassland_heather_km2),
    urban_km2 = as.numeric(urban_km2),
    suburban_km2 = as.numeric(suburban_km2),
    arable_km2 = as.numeric(arable_km2)
  ) %>%
  mutate(
    z_wetland   = scale(wetland_km2)[,1],
    z_saltwater = scale(saltwater_km2)[,1],
    z_freshwater = scale(freshwater_km2)[,1],
    z_woodland = scale(woodland_km2)[,1],
    z_grassland_heather = scale(grassland_heather_km2)[,1],
    z_urban = scale(urban_km2)[,1],
    z_suburban = scale(suburban_km2)[,1],
    z_arable = scale(arable_km2)[,1]
  )
z_woodland         <- landcover_grouped %>% dplyr::select(grid_id, z_woodland) %>% distinct() %>%arrange(grid_id) %>%pull(z_woodland)# create a vector
z_wetland          <- landcover_grouped %>% distinct(grid_id, z_wetland) %>% arrange(grid_id) %>% pull(z_wetland)
z_saltwater        <- landcover_grouped %>% distinct(grid_id, z_saltwater) %>% arrange(grid_id) %>% pull(z_saltwater)
z_freshwater       <- landcover_grouped %>% distinct(grid_id, z_freshwater) %>% arrange(grid_id) %>% pull(z_freshwater)
z_grassland_heather<- landcover_grouped %>% distinct(grid_id, z_grassland_heather) %>% arrange(grid_id) %>% pull(z_grassland_heather)
z_urban            <- landcover_grouped %>% distinct(grid_id, z_urban) %>% arrange(grid_id) %>% pull(z_urban)
z_suburban         <- landcover_grouped %>% distinct(grid_id, z_suburban) %>% arrange(grid_id) %>% pull(z_suburban)
z_arable           <- landcover_grouped %>% distinct(grid_id, z_arable) %>% arrange(grid_id) %>% pull(z_arable)

z_land <- cbind(
  z_wetland,
  z_freshwater,
  z_woodland,
  z_grassland_heather,
  z_arable,
  z_livestock,
  z_elevation
)

# z_land_df <- landcover_grouped %>%
#   distinct(grid_id, z_wetland, z_freshwater, z_woodland,
#            z_grassland_heather, z_arable, z_livestock, z_elevation) %>%
#   arrange(grid_id)


z_land_df <- as.data.frame(z_land)

z_land_df$grid_id <- landcover_grouped %>% 
  distinct(grid_id) %>% 
  arrange(grid_id) %>% 
  pull(grid_id)


write_csv(z_land_df, "01_data/covariates/z_land.csv")



# --- CITIZEN SCIENCE THINNING PROBABILITY ---
# including covariates that can affect where and 
# when people look and report mosquitoes — not the
# mosquito’s ecological preference
# ----------------------------------------------

## Points of Interest: local human activity metrics
poi <- st_read("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/POI/poi_uk.gpkg")
# write.csv(poi, "~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/POI/poi.csv")


library(dplyr)
main_categories <- poi %>%
  distinct(main_category)
main_categories <- data.frame(main_category = unique(poi$main_category))
# write.csv(main_categories, "~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/POI/poi_categories.csv")


# Define the relevant categories
relevant_categories <- c(
  # Highly Relevant
  "park", "beach", "trail", "nature_reserve", "hiking_trail", 
  "campground", "lake", "river", "forest", "botanical_garden",
  "national_park", "state_park", "playground", "dog_park", 
  "golf_course", "outdoor_gear",
  
  # Moderately Relevant
  "restaurant", "cafe", "pub", "bar", "ice_cream_shop", 
  "food_stand", "food_truck", "beer_garden",
  
  # Residential/Accommodation
  "hotel", "resort", "guest_house", "bed_and_breakfast", 
  "cottage", "accommodation", "apartments", "vacation_rental_agents"
)

# Subset the data
poi_subset <- poi[poi$main_category %in% relevant_categories, ]

poi_subset_df <- poi_subset %>%
  group_by(main_category) %>%
  summarize(count = n())



grid_poi <- grid_clipped
grid_poi <- st_transform(grid_poi, crs = st_crs(poi)) # must be same CRS
st_crs(grid_poi)

# Perform spatial join to count POIs in each grid cell
poi_counts <- st_join(poi, grid_poi, join = st_within) # can also use poi_subset instead of poi

# Count POIs per grid cell
grid_poi_counts <- poi_counts %>%
  st_drop_geometry() %>%  # Remove geometry for counting
  group_by(grid_id) %>%   
  summarise(poi_count = n())

# Join the counts back to the grid
grid_poi <- grid_poi %>%
  left_join(grid_poi_counts, by = "grid_id")

# Replace NA values with 0 (for grid cells with no POIs)
grid_poi$poi_count[is.na(grid_poi$poi_count)] <- 0
# group anything over 100??
summary(grid_poi$poi_count)
# Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
# 0.000    0.000    0.000    4.607    1.000 4409.000 

grid_poi$poi_count_grouped <- ifelse(grid_poi$poi_count > 100, 100, grid_poi$poi_count)

summary(grid_poi$poi_count_grouped)



# library(ggplot2)
# ggplot() +
#   geom_sf(data = grid_poi, aes(fill = poi_count)) +
#   scale_fill_viridis_c(name = "POI Count") +
#   labs(title = "POI Counts by Grid Cell") +
#   theme_minimal()

clean_data <- grid_poi %>%
  st_drop_geometry() 

write.csv(clean_data, "01_data/covariates/poi_grouped_area2km.csv", row.names = FALSE)

# grid_poi <- read_csv("01_data/covariates/poi_grouped_area2km.csv")

# ============================================================================
# BUILD FINAL DATA SETS
# ============================================================================

make_covariates <- function(
    landcover_grouped,
    livestock_df,
    elevation_df,
    poi_df,
    building_df,
    temp_df,
    rain_df,
    land_out = "01_data/covariates/z_land.csv",
    climate_out = "01_data/covariates/z_climate.csv"
) {
  library(dplyr)
  library(readr)
  
  message("---- Building STATIC land covariates ----")
  
  # --- 1. Base landcover static variables ---
  static_df <- landcover_grouped %>%
    distinct(
      grid_id,
      wetland_km2, freshwater_km2, woodland_km2,
      grassland_heather_km2, arable_km2
    ) %>%
    arrange(grid_id)
  
  # --- Merge livestock ---
  static_df <- static_df %>%
    left_join(livestock_df %>% distinct(grid_id, livestock_density), by = "grid_id")
  
  # --- Merge elevation ---
  static_df <- static_df %>%
    left_join(elevation_df %>% distinct(grid_id, elevation), by = "grid_id")
  
  # --- Merge POI density ---
  static_df <- static_df %>%
    left_join(poi_df %>% distinct(grid_id, poi_count_grouped), by = "grid_id")
  
  # --- Merge building density ---
  static_df <- static_df %>%
    left_join(building_df %>% distinct(grid_id, building_count), by = "grid_id")
  
  # ---- AUTOMATIC Z-SCORING FOR STATIC VARS ----
  static_numeric <- static_df %>% dplyr::select(-grid_id)
  
  static_scaled <- as.data.frame(scale(static_numeric))
  
  z_land_df <- bind_cols(
    grid_id = static_df$grid_id,
    static_scaled %>% 
      rename_with(~ paste0("z_", .x))
  )
  
  # ---- WRITE OUTPUT ----
  write_csv(z_land_df, land_out)
  message(paste("Wrote:", land_out))
  
  
  message("---- Building TIME-VARYING climate covariates ----")
  
  # --- 2. Combine climate variables ---
  climate_df <- temp_df %>%
    dplyr::select(grid_id_area2km, date, min_temp_21d_celsius) %>%
    left_join(
      rain_df %>% dplyr::select(grid_id_area2km, date, rainfall_28d),
      by = c("grid_id_area2km", "date")
    ) %>%
    arrange(grid_id_area2km, date)
  
  # ---- AUTOMATIC Z-SCORING (panel data) ----
  climate_scaled <- climate_df %>%
    group_by(grid_id_area2km) %>%       # z-scored within grid_id
    mutate(
      z_temp_min = scale(min_temp_21d_celsius)[,1],
      z_rain     = scale(rainfall_28d)[,1]
    ) %>%
    ungroup() %>%
    dplyr::select(grid_id_area2km, date, z_temp_min, z_rain)
  
  # ---- WRITE OUTPUT ----
  write_csv(climate_scaled, climate_out)
  message(paste("Wrote:", climate_out))
  
  invisible(
    list(
      z_land = z_land_df,
      z_climate = climate_scaled
    )
  )
}

make_covariates(
  landcover_grouped      = landcover_grouped,
  livestock_df           = grid_area2km_livestock,
  elevation_df           = grid_area2km_elevation2,
  poi_df                 = grid_poi,
  building_df            = building_data,
  temp_df                = temp_2km,
  rain_df                = rainfall_area2km
)



