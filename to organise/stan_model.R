# ==============================================================================
# STAN MARGINILIZATION OF ThE LATENT DISCRETE PARAMETER (N) USING A POISSON
# DISTRIBUTION INSTEAD OF A NEGATIVE BINOMIAL
# ==============================================================================
# ==============================================================================
# CMDSTAN WORKFLOW FOR MARGINALIZED POISSON MODEL
# ==============================================================================
# ==============================================================================
# 0. LIBRARY AND DATA LOADING
# ==============================================================================

setwd("/users/2601581k/pipiens-ISDM")

cat("Loading libraries...\n")

required_packages <- c(
  "readr", "writexl", "readxl", "tidyverse", "dplyr", "lubridate", "nimble",
  "sp", "Matrix", "pbapply"
)

# install.packages('spDataLarge', repos='https://nowosad.github.io/drat/', type='source')

install_missing <- function(packages) {
  installed <- packages %in% installed.packages()[,"Package"]
  if (any(!installed)) {
    install.packages(packages[!installed], repos = "https://cloud.r-project.org")
  }
}
install_missing(required_packages)

suppressPackageStartupMessages({
  library(readr)
  library(writexl)
  library(readxl)
  library(tidyverse)
  library(dplyr)
  library(lubridate)
  library(nimble)
  library(Matrix)
  library(sp)
  library(pbapply)
})


# ============================================================================
# LOAD COVARIATES (z-scored)
# ============================================================================
cat("Loading data...\n")
loaded_data <- load("data/loaded_data.RData")
z_land_data <- load("data/env_land_data.RData")
z_temp <- load("data/env_temp_data.RData")
z_rain_part1 <- load("data/env_rain_data_part1.RData")
z_rain_part2 <- load("data/env_rain_data_part2.RData")

load("data/loaded_data.RData")
load("data/env_land_data.RData")
load("data/env_temp_data.RData")
load("data/env_rain_data_part1.RData")
load("data/env_rain_data_part2.RData")
z_rain <- rbind(z_rain_part1, z_rain_part2)


remove(z_rain_part1)
remove(z_rain_part2)
gc()


# organise land covariates
z_land <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id_area2km) %>%
  dplyr::select(starts_with("z_")) %>%
  dplyr::select(-c(z_urban_km2, z_suburban_km2, z_poi_count_grouped, z_reports)) %>%
  as.matrix()

z_poi <- as.data.frame(z_land_data) %>%
  arrange(grid_id_area2km) %>%
  pull(z_poi_count_grouped)

z_reports <- as.data.frame(z_land_data) %>%
  arrange(grid_id_area2km) %>%
  pull(z_reports)


# ============================================================================
# CONVERT CLIMATE DATA TO MATRICES (grids × dates)
# ============================================================================

cat("\nConverting climate data to matrix format...\n")

# Get sorted unique dates for column ordering
# climate_dates_sorted <- sort(unique(z_climate_data$date))

# Create temperature matrix (rows = grids, columns = dates)
# z_temp <- z_climate_data %>%
#   dplyr::select(grid_id_area2km, date, z_temp_min) %>%
#   pivot_wider(
#     names_from = date,
#     values_from = z_temp_min
#   ) %>%
#   arrange(grid_id_area2km) %>%
#   dplyr::select(-grid_id_area2km) %>%
#   as.matrix()

# Create rainfall matrix (rows = grids, columns = dates)
# z_rain <- z_climate_data %>%
#   dplyr::select(grid_id_area2km, date, z_rain) %>%
#   pivot_wider(
#     names_from = date,
#     values_from = z_rain
#   ) %>%
#   arrange(grid_id_area2km) %>%
#   dplyr::select(-grid_id_area2km) %>%
#   as.matrix()

cat(sprintf("  Temperature matrix: %d grids × %d dates\n", nrow(z_temp), ncol(z_temp)))
cat(sprintf("  Rainfall matrix: %d grids × %d dates\n", nrow(z_rain), ncol(z_rain)))
cat(sprintf("  Temperature completeness: %.1f%%\n", 100 * mean(!is.na(z_temp))))
cat(sprintf("  Rainfall completeness: %.1f%%\n", 100 * mean(!is.na(z_rain))))


# Choose numerical constant for stability
# Larger values (1000-10000) help when probabilities are very small
CONSTANT <- 10000

cat("Presence-only records:", nrow(cs_presences), "\n")
cat("Using Bernoulli ones trick with CONSTANT =", CONSTANT, "\n\n")
