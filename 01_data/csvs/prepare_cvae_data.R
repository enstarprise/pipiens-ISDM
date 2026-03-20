
MODEL_TYPE <- "NNGP"   # "HSGP" | "NNGP" | "NNGP_centered" | "nonspatial"
GRID_SCOPE <- "all"        # "all"  | "observed"
QUADRATIC  <- "both"            # "none" | "rain" | "both"

stopifnot(MODEL_TYPE %in% c("HSGP", "NNGP", "NNGP_centered", "nonspatial"))
stopifnot(GRID_SCOPE %in% c("all", "observed"))
stopifnot(QUADRATIC  %in% c("none", "rain", "both"))
cat(sprintf("\n=== Spatial model: %s | Grid scope: %s | Quadratic: %s ===\n\n",
            MODEL_TYPE, GRID_SCOPE, QUADRATIC))

suppressPackageStartupMessages({
  library(readr)
  library(writexl)
  library(readxl)
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(lubridate)
  library(cmdstanr)
  library(Matrix)
  library(sp)
  library(sf)
  library(spdep)
  library(pbapply)
})

# ============================================================================
# 1. LOAD DATA
# ============================================================================
cat("Loading data...\n")

load("01_data/loaded_data_cov_2point5km.RData")

z_land <- as.data.frame(z_land_data_clean) %>%
  dplyr::arrange(grid_id) %>%
  dplyr::select(starts_with("z_")) %>%
  dplyr::select(-c(z_saltwater_km2, z_urban_km2, z_suburban_km2,
                   z_poi_count, z_reports, z_buildings_count,
                   z_poi_log, z_livestock_density, z_freshwater_km2_log)) %>%
  as.matrix()

z_poi <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>%
  pull(z_poi_log)

z_reports <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>%
  pull(z_reports)

all_grid_ids <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>%
  pull(grid_id)

cat(sprintf("Full data: %d grids x %d dates\n", nrow(z_temp_clean), ncol(z_temp_clean)))

# ============================================================================
# 2. GRID SCOPE: SELECT ROWS TO USE
# ============================================================================

if (GRID_SCOPE == "observed") {
  
  cat("\nSubsetting to observed grids only...\n")
  
  observed_grid_ids <- sort(unique(c(survey_df$grid_id, cs_df$grid_id)))
  n_grids_total     <- length(observed_grid_ids)
  grid_positions    <- match(observed_grid_ids, all_grid_ids)
  
  stopifnot(!any(is.na(grid_positions)))
  
  cat(sprintf("  %d -> %d grids (%.1f%% reduction)\n",
              length(all_grid_ids), n_grids_total,
              100 * (1 - n_grids_total / length(all_grid_ids))))
  
  # Grid ID -> sequential index mapping for this scope
  grid_id_to_idx <- setNames(seq_along(observed_grid_ids), observed_grid_ids)
  
} else {  # "all"
  
  cat("\nUsing all grids...\n")
  
  observed_grid_ids <- all_grid_ids
  n_grids_total     <- length(all_grid_ids)
  grid_positions    <- seq_along(all_grid_ids)   # identity
  
  grid_id_to_idx <- setNames(seq_along(all_grid_ids), all_grid_ids)
  
  cat(sprintf("  %d grids\n", n_grids_total))
}

# Subset all grid-indexed arrays to the selected scope
coords_scope    <- coords_2point5km[grid_positions, , drop = FALSE]
z_temp_scope    <- z_temp_clean[grid_positions, , drop = FALSE]
z_rain_scope    <- z_rain_clean[grid_positions, , drop = FALSE]
z_land_scope    <- z_land[grid_positions, , drop = FALSE]
z_poi_scope     <- z_poi[grid_positions]
z_reports_scope <- z_reports[grid_positions]
area_grid_scope <- area_grid[grid_positions]

# Map observation-level grid IDs to sequential indices within selected scope
survey_df_stan <- survey_df %>%
  mutate(grid_idx = grid_id_to_idx[as.character(grid_id)]) %>%
  filter(!is.na(Culex), !is.na(grid_idx))

cs_presences_stan <- cs_presences %>%
  mutate(grid_idx = grid_id_to_idx[as.character(grid_id)]) %>%
  filter(!is.na(grid_idx))

cat(sprintf("  Survey obs: %d | PO records: %d\n",
            nrow(survey_df_stan), nrow(cs_presences_stan)))

# ============================================================================
# 3. MISSING VALUE IMPUTATION
# ============================================================================
cat("\nHandling missing values...\n")

impute_missing_values <- function(x, method = "mean", verbose = TRUE) {
  if (!any(is.na(x))) {
    if (verbose) cat("  No missing values\n")
    return(x)
  }
  na_count <- sum(is.na(x))
  na_pct   <- 100 * na_count / length(x)
  if (method == "mean") {
    if (is.matrix(x)) {
      for (i in 1:nrow(x)) {
        row_na <- is.na(x[i, ])
        if (any(row_na) && !all(row_na)) x[i, row_na] <- mean(x[i, ], na.rm = TRUE)
      }
      for (j in 1:ncol(x)) {
        col_na <- is.na(x[, j])
        if (any(col_na) && !all(col_na)) x[col_na, j] <- mean(x[, j], na.rm = TRUE)
      }
      if (any(is.na(x))) x[is.na(x)] <- mean(x, na.rm = TRUE)
    } else {
      x[is.na(x)] <- mean(x, na.rm = TRUE)
    }
  } else if (method == "zero") {
    x[is.na(x)] <- 0
  }
  if (verbose) cat(sprintf("  Imputed %d values (%.1f%%) using %s\n", na_count, na_pct, method))
  return(x)
}

cat("Temperature:\n");    z_temp_scope    <- impute_missing_values(z_temp_scope)
cat("Rainfall:\n");       z_rain_scope    <- impute_missing_values(z_rain_scope)
cat("Land covariates:\n"); z_land_scope   <- impute_missing_values(z_land_scope)
cat("Obs-level RH:\n");   z_RH_clean      <- impute_missing_values(survey_df$z_RH)
cat("Obs-level rain:\n"); z_WS_rain_clean <- impute_missing_values(survey_df$z_WS_rain)

stopifnot(!any(is.na(z_temp_scope)), !any(is.na(z_rain_scope)),
          !any(is.na(z_land_scope)), !any(is.na(z_RH_clean)),
          !any(is.na(z_WS_rain_clean)))
cat("All missing values handled\n\n")

# ============================================================================
# 3B. QUADRATIC TERMS
# ============================================================================
# Quadratic matrices are built by squaring the already-standardised z-scores.
# They are SEPARATE from the linear matrices — both terms enter the model
# simultaneously. 


zeros_scope <- matrix(0, nrow = nrow(z_temp_scope), ncol = ncol(z_temp_scope))

z_temp_sq_scope <- if (QUADRATIC == "both") z_temp_scope^2 else zeros_scope
z_rain_sq_scope <- if (QUADRATIC %in% c("rain", "both")) z_rain_scope^2 else zeros_scope

cat(sprintf("Quadratic terms: %s\n\n", switch(QUADRATIC,
                                              "none" = "none  (z_temp_sq and z_rain_sq = zero matrices)",
                                              "rain" = "rain  (z_rain_sq = z_rain^2; z_temp_sq = 0)",
                                              "both" = "both  (z_temp_sq = z_temp^2, z_rain_sq = z_rain^2)"
)))






#######
load("01_data/loaded_data_cov_2point5km.RData")

grid <- read_csv("01_data/csvs/grid_clipped_2point5km_wkt.csv")
grid <- grid[,c("grid_id", "area")]

cvae_grid <- read_csv("01_data/csvs/priorCVAE_grid_coords.csv")

cvae_grid <- grid %>%
  inner_join(cvae_grid, grid, by = "grid_id")

write_csv(cvae_grid, "01_data/csvs/priorCVAE_grid_coords_area.csv")


# temp <- read_csv("01_data/csvs/tasmean_grid_21d_celsius.csv")


cs_df <- cs_df[,c("presence", "grid_id","grid_id_1",  "date", "date_1", "area")]
survey_df <- survey_df[, c("Culex",  "z_WS_rain","z_RH" , "grid_id", "area" ,
                           "grid_id_1","date_1", "date", "trap_type_idx")]

# fill NAs with 0 
survey_df$z_RH[is.na(survey_df$z_RH)] <- 0
survey_df$z_WS_rain[is.na(survey_df$z_WS_rain)] <- 0

write_csv(survey_df,"priorCVAE_survey_df.csv")
write_csv(cs_df,"priorCVAE_cs_df.csv")


z_land <- as.data.frame(z_land_data_clean[, c("grid_id","z_wetland_km2", "z_woodland_km2",
                                "z_freshwater_km2", "z_arable_km2",
                                "z_grassland_heather_km2", "z_livestock_density",
                                "z_elevation", "z_buildings_count")])

z_thin <- as.data.frame(z_land_data_clean[,c("grid_id","z_poi_count", "z_reports")])

write_csv(z_land, "priorCVAE_z_land.csv")
write_csv(z_thin, "priorCVAE_thinning.csv")


# SPATIOTEMPROAL DATAFRAMES: Z_TEMP AND Z_RAIN
library(tidyr)
library(dplyr)

# Create date character column
z_climate_data <- z_climate_data %>%
  mutate(date_char = as.character(date))

# Create wide dataframe for z_temp_min
z_temp_wide <- z_climate_data %>%
  select(grid_id, date_char, z_temp_min) %>%
  pivot_wider(
    id_cols = grid_id,
    names_from = date_char,
    values_from = z_temp_min
  )

# Fill NAs with 0 for all columns except grid_id
z_temp_wide <- z_temp_wide %>%
  mutate(across(-grid_id, ~ ifelse(is.na(.), 0, .)))

# Create wide dataframe for z_rain
z_rain_wide <- z_climate_data %>%
  select(grid_id, date_char, z_rain) %>%
  pivot_wider(
    id_cols = grid_id,
    names_from = date_char,
    values_from = z_rain
  )

# Fill NAs with 0 for all columns except grid_id
z_rain_wide <- z_rain_wide %>%
  mutate(across(-grid_id, ~ ifelse(is.na(.), 0, .)))

# Sort dates in chronological order (optional but recommended)
date_cols <- sort(as.Date(colnames(z_temp_wide)[-1]))  # exclude grid_id
date_cols_char <- as.character(date_cols)

# Reorder columns
z_temp_wide <- z_temp_wide %>%
  select(grid_id, all_of(date_cols_char))

z_rain_wide <- z_rain_wide %>%
  select(grid_id, all_of(date_cols_char))

# Check dimensions
dim(z_temp_wide)  # Should be (14642, 179) - grid_id + 178 dates
dim(z_rain_wide)  # Should be (14642, 179) - grid_id + 178 dates

# Write to CSV files
write.csv(z_temp_wide, "priorCVAE_z_temp_wide.csv", row.names = FALSE)
write.csv(z_rain_wide, "priorCVAE_z_rain_wide.csv", row.names = FALSE)

# Check first few rows and columns
print("z_temp_wide preview:")
print(z_temp_wide[1:5, 1:5])

print("z_rain_wide preview:")
print(z_rain_wide[1:5, 1:5])



z_temp <- as.data.frame(z_temp_clean)
z_rain <- as.data.frame(z_rain_clean)

write_csv(z_temp, "priorCVAE_z_temp.csv")
write_csv(z_rain, "priorCVAE_z_rain.csv")

