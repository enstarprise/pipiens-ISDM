rm(list = ls())
gc()

library(terra)
library(sf)
library(future.apply)
library(dplyr)
library(exactextractr)
library(tidyr)

# --- scotland boundary & grid---
scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp")
scotland <- st_transform(scotland,crs = 27700)

grid_clipped <- st_read("grids/grid_clipped_1km.gpkg")

landcover <- rast("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/data/envData/landCover/UK_CEH/data/7727ce7d-531e-4d77-b756-5cc59ff016bd/ukregion-scotland.tif")

# grid at the crs of the landcover file
grid_for_extraction <- st_transform(grid_clipped, crs(landcover))
st_crs(landcover) == st_crs(grid_for_extraction)

get_landcover_area_per_class <- function(
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
  
  # Ensure grid has a grid_id
  if (!"grid_id" %in% names(grid_sf)) {
    grid_sf$grid_id <- seq_len(nrow(grid_sf))
  }
  
  # Load raster to get CRS and resolution
  base_raster <- terra::rast(landcover_raster_path)
  
  # Reproject grid if needed
  if (!identical(crs(grid_sf), crs(base_raster))) {
    message("Reprojecting grid to match raster CRS...")
    grid_sf <- st_transform(grid_sf, crs(base_raster))
  }
  
  pixel_area_km2 <- prod(terra::res(base_raster)) / 1e6
  
  # Parallel setup
  if (num_cores > 1) {
    future::plan(future::multisession, workers = num_cores)
  } else {
    future::plan(future::sequential)
  }
  
  # Chunk the grid
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
          
          # print(paste("Cell", j, "extracted rows:", if (is.null(extracted)) 0 else nrow(extracted)))
          # print(names(extracted))
          
          # # If extracted is null or zero rows, skip this polygon early
          # if (is.null(extracted) || nrow(extracted) == 0) return(NULL)
          
          # Safeguard against NULL or missing 'value'
          if (is.null(extracted) || !"Class" %in% names(extracted) || nrow(extracted) == 0) {
            return(NULL)
          }
          
          # Inside the loop
          extracted %>%
            filter(!is.na(Class)) %>%
            group_by(Class) %>%
            summarise(area_km2 = sum(coverage_fraction * pixel_area_km2, na.rm = TRUE), .groups = "drop") %>%
            mutate(grid_id = grid_sf$grid_id[j])
        }, error = function(e) {
          message(sprintf("Error in cell %d: %s", j, e$message))
          NULL
        })
      },
      future.seed = TRUE
    )
    
    # Save this chunk
    chunk_df <- bind_rows(chunk_results)
    
    saveRDS(
      chunk_df,
      file.path(output_dir, sprintf("chunk_%03d.rds", i))
    )
    
    rm(chunk_results)
    gc()
  }
  
  # Combine results
  result_files <- list.files(output_dir, pattern = "^chunk_.*\\.rds$", full.names = TRUE)
  all_results <- bind_rows(lapply(result_files, readRDS))
  
  # Cleanup temp files
  unlink(output_dir, recursive = TRUE)
  print(head(all_results))
  print(names(all_results))
  
  # Pivot to wide format
  wide_result <- all_results %>%
    mutate(class = paste0("class_", .data[["Class"]], "_km2")) %>%
    select(grid_id, class, area_km2) %>%
    tidyr::pivot_wider(
      names_from = class,
      values_from = area_km2,
      values_fill = 0
    )
  
  
  return(wide_result)
}

### --- EXTRACTION --- ###

# Save raster to disk before calling function
terra::writeRaster(landcover, "landcover.tif", overwrite = TRUE)

wide_class_areas <- get_landcover_area_per_class(
  grid_sf = grid_for_extraction,
  landcover_raster_path = "landcover.tif",
  num_cores = 10,
  chunk_size = 2000
)

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

write.csv(wide_class_areas, "csvs/landcover_1km.csv", row.names = FALSE)



