library(dplyr)

load("01_data/processedCovariates/1km/stan_data_init_1km_full.RData")



init_fun <- function() {
  list(
    # --------------------
    # Abundance model
    # --------------------
    beta0      = rnorm(1, 0, 0.2),
    beta_temp  = rnorm(1, 0, 0.2),
    beta_rain  = rnorm(1, 0, 0.2),
    beta_land  = rnorm(n_land_covs, 0, 0.2),
    
    # --------------------
    # Detection model
    # --------------------
    alpha0         = qlogis(0.1),
    alpha_RH       = rnorm(1, 0, 0.2),
    alpha_WS_rain  = rnorm(1, 0, 0.2),
    alpha_trap_raw = rnorm(4, 0, 0.1),
    log_sigma_trap = log(0.5),
    
    # --------------------
    # Thinning model
    # --------------------
    delta0        = rnorm(1, -1, 0.2),
    delta_poi     = rnorm(1, 0, 0.2),
    delta_reports = rnorm(1, 0, 0.2),
    
    # --------------------
    # Dispersion
    # --------------------
    phi = abs(rnorm(1, 2, 0.5)) + 0.5
  )
}

# Assuming you have your original data loaded
# original_stan_data <- your full dataset
subset_to_5k_pixels <- function(original_data) {
  
  # Step 1: Identify pixels that MUST be included (have survey or PO data)
  pixels_with_survey <- unique(original_data$survey_grid_idx)
  pixels_with_po <- unique(original_data$po_grid_idx)
  
  # Union of pixels with data
  required_pixels <- unique(c(pixels_with_survey, pixels_with_po))
  
  n_required <- length(required_pixels)
  
  cat(sprintf("Pixels with survey data: %d\n", length(pixels_with_survey)))
  cat(sprintf("Pixels with PO data: %d\n", length(pixels_with_po)))
  cat(sprintf("Total required pixels: %d\n", n_required))
  
  # Step 2: Check if we need additional pixels
  if (n_required > 5000) {
    stop(sprintf("Required pixels (%d) exceeds 5000! Cannot subset.", n_required))
  }
  
  n_additional_needed <- 5000 - n_required
  cat(sprintf("Additional pixels needed: %d\n", n_additional_needed))
  
  # Step 3: Sample additional pixels from remaining pixels
  all_pixels <- 1:original_data$n_grids_total
  remaining_pixels <- setdiff(all_pixels, required_pixels)
  
  if (n_additional_needed > 0) {
    additional_pixels <- sample(remaining_pixels, n_additional_needed)
    selected_pixels <- c(required_pixels, additional_pixels)
  } else {
    selected_pixels <- required_pixels
  }
  
  # Sort for easier indexing
  selected_pixels <- sort(selected_pixels)
  
  cat(sprintf("Total selected pixels: %d\n", length(selected_pixels)))
  
  # Step 4: Create mapping from old grid ID to new grid ID (1:length(selected_pixels))
  old_to_new <- setNames(1:length(selected_pixels), selected_pixels)
  
  # Step 5: Remap survey and PO grid indices to the new subset
  new_survey_grid_idx <- old_to_new[as.character(original_data$survey_grid_idx)]
  new_po_grid_idx <- old_to_new[as.character(original_data$po_grid_idx)]
  
  # Step 6: Define obs_grid_idx as all grids that have survey observations
  new_obs_grid_idx <- sort(unique(new_survey_grid_idx))
  
  # Step 7: site_to_grid indexes each survey observation within obs_grid_idx
  new_site_to_grid <- match(new_survey_grid_idx, new_obs_grid_idx)
  
  # Step 8: Create the subset stan data
  stan_data_5k <- list(
    # Dimensions
    n_grids_total = length(selected_pixels),
    n_grids_obs = length(new_obs_grid_idx),
    n_dates = original_data$n_dates,
    n_obs_y = original_data$n_obs_y,
    n_obs_po = original_data$n_obs_po,
    n_land_covs = original_data$n_land_covs,
    
    # Indices - REMAPPED
    obs_grid_idx = new_obs_grid_idx,
    site_to_grid = new_site_to_grid,
    date_y = original_data$date_y,
    trap_type = original_data$trap_type,
    po_grid_idx = new_po_grid_idx,
    survey_grid_idx = new_survey_grid_idx,
    date_po = original_data$date_po,
    
    # Observations (unchanged)
    y = original_data$y,
    ones = original_data$ones,
    
    # Covariates - SUBSET to selected pixels
    z_temp = original_data$z_temp[selected_pixels, , drop = FALSE],
    z_rain = original_data$z_rain[selected_pixels, , drop = FALSE],
    z_land = original_data$z_land[selected_pixels, , drop = FALSE],
    z_poi = original_data$z_poi[selected_pixels],
    z_reports = original_data$z_reports[selected_pixels],
    
    # Site-level covariates (unchanged - per observation)
    z_RH = original_data$z_RH,
    z_WS_rain = original_data$z_WS_rain,
    
    # Other - SUBSET
    area_grid = original_data$area_grid[selected_pixels],
    CONSTANT = original_data$CONSTANT,
    N_multiplier = original_data$N_multiplier
  )
  
  # Step 9: Validation checks
  cat("\n=== Validation ===\n")
  cat(sprintf("All site_to_grid indices valid: %s\n", 
              all(!is.na(stan_data_5k$site_to_grid) & 
                    stan_data_5k$site_to_grid >= 1 & 
                    stan_data_5k$site_to_grid <= stan_data_5k$n_grids_obs)))
  
  cat(sprintf("All survey_grid_idx indices valid: %s\n",
              all(!is.na(stan_data_5k$survey_grid_idx) & 
                    stan_data_5k$survey_grid_idx >= 1 & 
                    stan_data_5k$survey_grid_idx <= stan_data_5k$n_grids_total)))
  
  cat(sprintf("All po_grid_idx indices valid: %s\n",
              all(!is.na(stan_data_5k$po_grid_idx) & 
                    stan_data_5k$po_grid_idx >= 1 & 
                    stan_data_5k$po_grid_idx <= stan_data_5k$n_grids_total)))
  
  cat(sprintf("site_to_grid correctly maps into obs_grid_idx: %s\n",
              all(stan_data_5k$obs_grid_idx[stan_data_5k$site_to_grid] == stan_data_5k$survey_grid_idx)))
  
  # Return both the data and the mapping for reference
  return(list(
    stan_data = stan_data_5k,
    selected_pixels = selected_pixels,
    old_to_new_mapping = old_to_new
  ))
}



# Usage
result <- subset_to_5k_pixels(stan_data)

# Print summary
cat("\n=== Summary ===\n")
cat(sprintf("Original pixels: %d\n", stan_data$n_grids_total))
cat(sprintf("Subset pixels: %d\n", result$stan_data$n_grids_total))
cat(sprintf("Survey observations: %d\n", result$stan_data$n_obs_y))
cat(sprintf("PO observations: %d\n", result$stan_data$n_obs_po))

# Additional debugging
cat("\n=== Sample Mappings ===\n")
cat("Original survey_grid_idx (first 10):", head(stan_data$survey_grid_idx, 10), "\n")
cat("New survey_grid_idx (first 10):", head(result$stan_data$survey_grid_idx, 10), "\n")
cat("Original po_grid_idx (first 10):", head(stan_data$po_grid_idx, 10), "\n")
cat("New po_grid_idx (first 10):", head(result$stan_data$po_grid_idx, 10), "\n")

# Save 
saveRDS(result$stan_data, "stan_data_5k.rds")
saveRDS(result$selected_pixels, "selected_pixels_5k.rds")
saveRDS(result$old_to_new_mapping, "pixel_mapping_5k.rds")

subset_stan_data <- result$stan_data
selected_pixels <- result$selected_pixels
old_to_new_mapping <- result$old_to_new_mapping


save(subset_stan_data, selected_pixels, old_to_new_mapping, result, init_fun,
     n_dates, n_land_covs, observed_grid_ids,
     file = "02_model/stan/01_5kpixel_test/5k_pixelSubset.RData")

