
# =============================================================================
# COMPUTE HSGP BOUNDARIES FROM SHAPEFILE
# Reduce L_x, L_y for HSGP without removing any grid cells from model
# =============================================================================

library(sf)

# =============================================================================
# MAIN FUNCTION: COMPUTE L FROM BOUNDARY SHAPEFILE
# =============================================================================

compute_L_from_boundary <- function(grid_coords_centered, boundary_sf, 
                                    grid_sf = NULL, c_factor = 1.5) {
  #' Compute HSGP L values based on boundary shapefile extent
  #' Does NOT filter grid cells - just computes tighter boundaries for HSGP
  #' 
  #' @param grid_coords_centered: centered grid coordinates (matrix, n×2)
  #' @param boundary_sf: sf object with study area boundary
  #' @param grid_sf: optional, sf object of grid (for CRS if needed)
  #' @param c_factor: boundary extension factor (typically 1.2-1.5)
  #' @return: list with L_x, L_y and diagnostic information
  
  # Get bounding box of boundary shapefile
  bbox <- st_bbox(boundary_sf)
  
  # Extract coordinates
  # bbox format: xmin, ymin, xmax, ymax
  x_min <- bbox["xmin"]
  x_max <- bbox["xmax"]
  y_min <- bbox["ymin"]
  y_max <- bbox["ymax"]
  
  # Convert to km if needed (check if values > 10000, suggesting meters)
  if (max(abs(c(x_min, x_max, y_min, y_max))) > 10000) {
    cat("Converting boundary coordinates from meters to kilometers...\n")
    x_min <- x_min / 1000
    x_max <- x_max / 1000
    y_min <- y_min / 1000
    y_max <- y_max / 1000
  }
  
  # Find the center of the boundary
  x_center_boundary <- (x_min + x_max) / 2
  y_center_boundary <- (y_min + y_max) / 2
  
  # Find the center of the grid
  x_center_grid <- mean(range(grid_coords_centered[,1]))
  y_center_grid <- mean(range(grid_coords_centered[,2]))
  
  cat("\n=== BOUNDARY SHAPEFILE ANALYSIS ===\n")
  cat("Boundary extent:\n")
  cat("  X: [", round(x_min, 2), ",", round(x_max, 2), "] km\n")
  cat("  Y: [", round(y_min, 2), ",", round(y_max, 2), "] km\n")
  cat("  Width:", round(x_max - x_min, 2), "km\n")
  cat("  Height:", round(y_max - y_min, 2), "km\n")
  
  cat("\nCenters:\n")
  cat("  Boundary center: (", round(x_center_boundary, 2), ",", 
      round(y_center_boundary, 2), ")\n")
  cat("  Grid center: (", round(x_center_grid, 2), ",", 
      round(y_center_grid, 2), ")\n")
  
  # Check if centers align (they should if grid_coords are properly centered)
  center_offset <- sqrt((x_center_boundary - x_center_grid)^2 + 
                          (y_center_boundary - y_center_grid)^2)
  
  if (center_offset > 10) {
    cat("  WARNING: Centers don't align! Offset =", round(center_offset, 2), "km\n")
    cat("  This suggests grid_coords may not be centered properly.\n")
  } else {
    cat("  Centers align well (offset =", round(center_offset, 2), "km)\n")
  }
  
  # Compute L based on boundary extent
  # L should be the maximum absolute coordinate from center
  # For a boundary centered at origin, this is max(abs(x_min), abs(x_max))
  
  # Re-center boundary coordinates to match grid center (0, 0)
  x_min_centered <- x_min - x_center_boundary
  x_max_centered <- x_max - x_center_boundary
  y_min_centered <- y_min - y_center_boundary
  y_max_centered <- y_max - y_center_boundary
  
  # L is the max absolute coordinate, then extended by c_factor
  L_x <- c_factor * max(abs(x_min_centered), abs(x_max_centered))
  L_y <- c_factor * max(abs(y_min_centered), abs(y_max_centered))
  
  # Also compute what L would be for the full grid (for comparison)
  L_x_full_grid <- c_factor * max(abs(range(grid_coords_centered[,1])))
  L_y_full_grid <- c_factor * max(abs(range(grid_coords_centered[,2])))
  
  cat("\n=== HSGP DOMAIN PARAMETERS ===\n")
  cat("Boundary factor (c):", c_factor, "\n\n")
  
  cat("ORIGINAL (full grid):\n")
  cat("  L_x:", round(L_x_full_grid, 2), "km\n")
  cat("  L_y:", round(L_y_full_grid, 2), "km\n")
  cat("  L (geometric mean):", round(sqrt(L_x_full_grid * L_y_full_grid), 2), "km\n")
  
  cat("\nNEW (boundary-based):\n")
  cat("  L_x:", round(L_x, 2), "km\n")
  cat("  L_y:", round(L_y, 2), "km\n")
  cat("  L (geometric mean):", round(sqrt(L_x * L_y), 2), "km\n")
  
  cat("\nREDUCTION:\n")
  cat("  L_x:", round(100 * (1 - L_x/L_x_full_grid), 1), "%\n")
  cat("  L_y:", round(100 * (1 - L_y/L_y_full_grid), 1), "%\n")
  
  # Check if any grid cells fall outside new boundaries
  outside_x <- sum(abs(grid_coords_centered[,1]) > L_x)
  outside_y <- sum(abs(grid_coords_centered[,2]) > L_y)
  outside_both <- sum(abs(grid_coords_centered[,1]) > L_x | 
                        abs(grid_coords_centered[,2]) > L_y)
  
  cat("\n=== GRID CELLS OUTSIDE NEW BOUNDARIES ===\n")
  cat("Cells beyond L_x:", outside_x, "\n")
  cat("Cells beyond L_y:", outside_y, "\n")
  cat("Cells beyond either:", outside_both, 
      "(", round(100 * outside_both / nrow(grid_coords_centered), 1), "%)\n")
  
  if (outside_both > 0) {
    cat("\nNOTE: These cells will remain in the model but will have\n")
    cat("      less accurate HSGP approximation (edge effects).\n")
    cat("      This is acceptable if they're outlying islands.\n")
  }
  
  list(
    L_x = L_x,
    L_y = L_y,
    L_x_full_grid = L_x_full_grid,
    L_y_full_grid = L_y_full_grid,
    boundary_extent_x = x_max - x_min,
    boundary_extent_y = y_max - y_min,
    n_cells_outside = outside_both,
    pct_cells_outside = 100 * outside_both / nrow(grid_coords_centered)
  )
}


# =============================================================================
# ESTIMATE M_SQRT WITH NEW L VALUES
# =============================================================================

estimate_M_sqrt_for_new_L <- function(L_x, L_y, rho_target = 4.0,
                                      kernel = "matern32", c_factor = 1.5) {
  #' Estimate required M_sqrt for new L values
  #' (Same as before, but separated for clarity)
  
  # Coefficients from Riutort-Mayol et al. (2023)
  coeffs <- list(
    matern32 = 3.42,
    matern52 = 2.65,
    squared_exp = 1.75,
    matern12 = 3.42  # Use matern32 as conservative estimate
  )
  
  coeff <- coeffs[[kernel]]
  
  # Compute minimum m per dimension
  m_x <- coeff * c_factor * (L_x / rho_target)
  m_y <- coeff * c_factor * (L_y / rho_target)
  
  M_sqrt_needed <- ceiling(max(m_x, m_y))
  
  # Practical values
  practical <- c(50, 64, 100, 128, 150, 200, 256, 300, 400, 500, 512, 
                 600, 700, 800, 900, 1000, 1024)
  M_sqrt_suggested <- practical[which.min(abs(practical - M_sqrt_needed))]
  
  # Minimum resolvable lengthscale with suggested M
  rho_min <- (coeff * c_factor * max(L_x, L_y)) / M_sqrt_suggested
  
  cat("\n=== M_SQRT ESTIMATION ===\n")
  cat("Target lengthscale (ρ):", rho_target, "km\n")
  cat("Kernel:", kernel, "(coefficient =", coeff, ")\n\n")
  
  cat("Required m per dimension:\n")
  cat("  m_x =", round(m_x, 1), "\n")
  cat("  m_y =", round(m_y, 1), "\n")
  cat("  M_sqrt needed:", M_sqrt_needed, "\n")
  cat("  M_sqrt suggested:", M_sqrt_suggested, "\n")
  cat("  Total M:", M_sqrt_suggested^2, "basis functions\n\n")
  
  cat("Diagnostic check:\n")
  cat("  With M_sqrt =", M_sqrt_suggested, "\n")
  cat("  Minimum resolvable ρ:", round(rho_min, 2), "km\n")
  
  if (rho_target >= rho_min) {
    cat("  ✓ Can accurately approximate ρ =", rho_target, "km\n")
  } else {
    cat("  ✗ Cannot accurately approximate ρ =", rho_target, "km\n")
    cat("    Need M_sqrt ≥", M_sqrt_needed, "\n")
  }
  
  list(
    M_sqrt_needed = M_sqrt_needed,
    M_sqrt_suggested = M_sqrt_suggested,
    M = M_sqrt_suggested^2,
    rho_min = rho_min,
    m_x = m_x,
    m_y = m_y
  )
}


# =============================================================================
# COMPLETE WORKFLOW
# =============================================================================

hsgp_boundaries_from_shapefile <- function(grid_coords_centered, 
                                           boundary_sf,
                                           c_factor = 1.5,
                                           rho_target = 4.0,
                                           kernel = "matern32") {
  #' Complete workflow: compute L and M_sqrt from boundary shapefile
  #' 
  #' @param grid_coords_centered: centered grid coordinates (matrix)
  #' @param boundary_sf: boundary shapefile
  #' @param c_factor: boundary extension factor
  #' @param rho_target: expected/target lengthscale
  #' @param kernel: kernel type
  #' @return: list with L values and M_sqrt recommendation
  
  cat("=================================================================\n")
  cat("HSGP DOMAIN BOUNDARY COMPUTATION FROM SHAPEFILE\n")
  cat("=================================================================\n")
  
  # Step 1: Compute L from boundary
  L_result <- compute_L_from_boundary(grid_coords_centered, boundary_sf, 
                                      c_factor = c_factor)
  
  # Step 2: Estimate M_sqrt
  M_result <- estimate_M_sqrt_for_new_L(L_result$L_x, L_result$L_y,
                                        rho_target, kernel, c_factor)
  
  cat("\n=================================================================\n")
  cat("SUMMARY\n")
  cat("=================================================================\n")
  cat("Update your Stan data with:\n")
  cat("  stan_data$L_x <-", round(L_result$L_x, 2), "\n")
  cat("  stan_data$L_y <-", round(L_result$L_y, 2), "\n")
  cat("  stan_data$M_sqrt <-", M_result$M_sqrt_suggested, "\n")
  cat("\nGrid cells remain unchanged.\n")
  cat("All data structures remain unchanged.\n")
  cat("Only the HSGP approximation domain is reduced.\n")
  cat("=================================================================\n")
  
  list(
    L_x = L_result$L_x,
    L_y = L_result$L_y,
    M_sqrt = M_result$M_sqrt_suggested,
    L_result = L_result,
    M_result = M_result
  )
}


# =============================================================================
# USAGE EXAMPLES
# =============================================================================

result <- hsgp_boundaries_from_shapefile(grid_coords_centered = grid_coords_centered,
      boundary_sf = st_read('01_data/grids/main_landmass.shp'), c_factor = 1.5,
      rho_target = 4.0, kernel = 'matern32')



cat("\n=== Script loaded successfully ===\n\n")
cat("USAGE:\n\n")
cat("# Simple usage:\n")
cat("result <- hsgp_boundaries_from_shapefile(\n")
cat("  grid_coords_centered = your_centered_coords,\n")
cat("  boundary_sf = st_read('mainland_scotland.shp'),\n")
cat("  c_factor = 1.5,\n")
cat("  rho_target = 4.0,\n")
cat("  kernel = 'matern32'\n")
cat(")\n\n")
cat("# Then update Stan data:\n")
cat("stan_data$L_x <- result$L_x\n")
cat("stan_data$L_y <- result$L_y\n")
cat("stan_data$M_sqrt <- result$M_sqrt\n\n")
cat("# That's it! No other changes needed.\n\n")




# =============================================================================
# TRIM HSGP DOMAIN TO MAIN DATA EXTENT
# Exclude outlying islands to reduce computational burden
# =============================================================================

library(sf)
library(dplyr)
library(dbscan)

# Assuming you have:
# - grid_sf: your spatial grid (sf object)
# - grid_coords_centered: centered coordinates matrix



# =============================================================================
# METHOD 1: USE QUANTILES (Simple and Fast)
# =============================================================================
# Exclude extreme 1-2% of points in each dimension

trim_domain_quantiles <- function(grid_coords_centered, prob = 0.01) {
  # Keep central 98% (exclude 1% on each tail)
  x_limits <- quantile(grid_coords_centered[,1], c(prob, 1-prob))
  y_limits <- quantile(grid_coords_centered[,2], c(prob, 1-prob))
  
  # Identify cells within limits
  keep_idx <- grid_coords_centered[,1] >= x_limits[1] & 
    grid_coords_centered[,1] <= x_limits[2] &
    grid_coords_centered[,2] >= y_limits[1] & 
    grid_coords_centered[,2] <= y_limits[2]
  
  # Get trimmed extents
  x_range <- range(grid_coords_centered[keep_idx, 1])
  y_range <- range(grid_coords_centered[keep_idx, 2])
  
  list(
    keep_idx = keep_idx,
    x_range = x_range,
    y_range = y_range,
    extent_x = diff(x_range),
    extent_y = diff(y_range),
    n_removed = sum(!keep_idx),
    pct_removed = 100 * mean(!keep_idx)
  )
}

# Example usage:
result <- trim_domain_quantiles(grid_coords_centered, prob = 0.01)
cat("Removed", result$n_removed, "cells (", round(result$pct_removed, 2), "%)\n")
# cat("New extent X:", result$extent_x, "km\n")
# cat("New extent Y:", result$extent_y, "km\n")


# =============================================================================
# METHOD 2: DENSITY-BASED TRIMMING (Better for Islands)
# =============================================================================
# Remove isolated cells using spatial clustering

trim_domain_density <- function(grid_sf, min_cluster_size = 50) {
  library(dbscan)
  
  # Get centroids
  centroids <- st_coordinates(st_centroid(grid_sf))
  
  # Cluster using DBSCAN (finds dense regions)
  # eps = max distance for neighbors (in km)
  # minPts = minimum points to form cluster
  clusters <- dbscan(centroids, eps = 50, minPts = 5)
  
  # Find largest cluster (mainland)
  cluster_sizes <- table(clusters$cluster)
  # Cluster 0 = noise/outliers, exclude it
  cluster_sizes <- cluster_sizes[names(cluster_sizes) != "0"]
  main_cluster <- as.numeric(names(which.max(cluster_sizes)))
  
  # Keep only main cluster
  keep_idx <- clusters$cluster == main_cluster
  
  # Also keep any cluster larger than min_cluster_size
  for (clust in names(cluster_sizes)) {
    if (cluster_sizes[clust] >= min_cluster_size) {
      keep_idx <- keep_idx | (clusters$cluster == as.numeric(clust))
    }
  }
  
  # Get trimmed coordinates
  coords_trimmed <- centroids[keep_idx, ]
  x_range <- range(coords_trimmed[, 1])
  y_range <- range(coords_trimmed[, 2])
  
  list(
    keep_idx = keep_idx,
    x_range = x_range,
    y_range = y_range,
    extent_x = diff(x_range),
    extent_y = diff(y_range),
    n_removed = sum(!keep_idx),
    pct_removed = 100 * mean(!keep_idx),
    n_clusters = length(unique(clusters$cluster[clusters$cluster != 0]))
  )
}

# Example usage:
result <- trim_domain_density(grid, min_cluster_size = 50)
cat("Found", result$n_clusters, "main cluster(s)\n")
cat("Removed", result$n_removed, "cells (", round(result$pct_removed, 2), "%)\n")


# =============================================================================
# METHOD 3: CONVEX HULL (Most Conservative)
# =============================================================================
# Keep everything within convex hull of actual observations

trim_domain_observations <- function(grid_sf, obs_locations_sf) {
  # Create convex hull around observations
  obs_union <- st_union(obs_locations_sf)
  hull <- st_convex_hull(obs_union)
  
  # Add buffer (e.g., 50 km) to include areas near observations
  hull_buffered <- st_buffer(hull, dist = 50000)  # 50 km in meters
  
  # Check which grid cells intersect the hull
  keep_idx <- st_intersects(grid_sf, hull_buffered, sparse = FALSE)[,1]
  
  # Get trimmed coordinates
  grid_trimmed <- grid_sf[keep_idx, ]
  centroids <- st_coordinates(st_centroid(grid_trimmed))
  
  x_range <- range(centroids[, 1])
  y_range <- range(centroids[, 2])
  
  list(
    keep_idx = keep_idx,
    x_range = x_range,
    y_range = y_range,
    extent_x = diff(x_range),
    extent_y = diff(y_range),
    n_removed = sum(!keep_idx),
    pct_removed = 100 * mean(!keep_idx)
  )
}

# Example usage (if you have observation locations):
# result <- trim_domain_observations(grid_sf, survey_sites_sf)


# =============================================================================
# METHOD 4: MANUAL INSPECTION AND TRIMMING
# =============================================================================

manual_trim_islands <- function(grid_coords_centered, grid_sf) {
  # Plot to identify outliers
  plot(grid_coords_centered[,1], grid_coords_centered[,2], 
       pch = 16, cex = 0.3,
       main = "Grid Cells - Identify Islands to Remove",
       xlab = "X (km)", ylab = "Y (km)")
  
  # Print extreme cells
  cat("\nMost extreme cells:\n")
  cat("\nFarthest WEST:\n")
  west_idx <- which.min(grid_coords_centered[,1])
  print(grid_coords_centered[west_idx,])
  
  cat("\nFarthest EAST:\n")
  east_idx <- which.max(grid_coords_centered[,1])
  print(grid_coords_centered[east_idx,])
  
  cat("\nFarthest SOUTH:\n")
  south_idx <- which.min(grid_coords_centered[,2])
  print(grid_coords_centered[south_idx,])
  
  cat("\nFarthest NORTH:\n")
  north_idx <- which.max(grid_coords_centered[,2])
  print(grid_coords_centered[north_idx,])
  
  # Identify cells beyond certain thresholds
  cat("\n\nCells distribution by quadrant:\n")
  cat("X > 400 km:", sum(grid_coords_centered[,1] > 400), "cells\n")
  cat("X < -400 km:", sum(grid_coords_centered[,1] < -400), "cells\n")
  cat("Y > 400 km:", sum(grid_coords_centered[,2] > 400), "cells\n")
  cat("Y < -200 km:", sum(grid_coords_centered[,2] < -200), "cells\n")
  
  # You can manually set boundaries based on plot
  # For example, if you see islands beyond certain limits:
  # x_min <- -400  # Exclude Hebrides if too far
  # x_max <- 300   # Exclude Shetland if too far  
  # y_min <- -200  # Southern boundary
  # y_max <- 450   # Northern boundary (exclude Shetland)
  
  # Return structure for you to fill in
  list(
    message = "Inspect plot and set manual boundaries"
  )
}

# Example usage:
manual_trim_islands(grid_coords_centered, grid)
# Then define boundaries manually


# =============================================================================
# RECOMMENDED WORKFLOW
# =============================================================================

recommended_workflow <- function(grid_coords_centered, grid_sf) {
  cat("=== HSGP DOMAIN TRIMMING WORKFLOW ===\n\n")
  
  # Step 1: Check current extent
  cat("CURRENT DOMAIN:\n")
  cat("X range:", range(grid_coords_centered[,1]), "km\n")
  cat("Y range:", range(grid_coords_centered[,2]), "km\n")
  cat("Extent X:", diff(range(grid_coords_centered[,1])), "km\n")
  cat("Extent Y:", diff(range(grid_coords_centered[,2])), "km\n\n")
  
  # Step 2: Try quantile trimming (removes outliers)
  cat("OPTION 1: Quantile trimming (exclude extreme 2%):\n")
  result1 <- trim_domain_quantiles(grid_coords_centered, prob = 0.01)
  cat("  New extent X:", result1$extent_x, "km\n")
  cat("  New extent Y:", result1$extent_y, "km\n")
  cat("  Cells removed:", result1$n_removed, "\n\n")
  
  # Step 3: Visualize
  cat("STEP 3: Visualize to identify islands:\n")
  plot(grid_coords_centered[,1], grid_coords_centered[,2], 
       pch = 16, cex = 0.3, col = "gray",
       main = "All Cells (gray) vs Trimmed (red)",
       xlab = "X (km)", ylab = "Y (km)")
  points(grid_coords_centered[result1$keep_idx, 1],
         grid_coords_centered[result1$keep_idx, 2],
         pch = 16, cex = 0.3, col = "red")
  
  cat("\nInspect the plot to decide on final trimming strategy.\n")
  
  return(result1)
}


# =============================================================================
# COMPUTE NEW L VALUES AFTER TRIMMING
# =============================================================================

compute_new_L <- function(trim_result, c_factor = 1.5) {
  # L should be c × max absolute coordinate
  L_x <- c_factor * max(abs(trim_result$x_range))
  L_y <- c_factor * max(abs(trim_result$y_range))
  
  cat("\n=== NEW HSGP PARAMETERS ===\n")
  cat("L_x:", L_x, "km\n")
  cat("L_y:", L_y, "km\n")
  
  # Estimate M_sqrt needed for rho ~ 4 km
  rho_target <- 4.0
  m_x <- 3.42 * c_factor * (L_x / rho_target)
  m_y <- 3.42 * c_factor * (L_y / rho_target)
  M_sqrt_needed <- ceiling(max(m_x, m_y))
  
  cat("\nFor ρ ≈", rho_target, "km (Matérn 3/2):\n")
  cat("  M_sqrt needed:", M_sqrt_needed, "\n")
  cat("  Total M:", M_sqrt_needed^2, "basis functions\n")
  
  # Suggest practical M_sqrt
  practical_M <- c(256, 400, 500, 600, 800, 1000)
  suggested <- practical_M[which.min(abs(practical_M - M_sqrt_needed))]
  
  cat("\nSuggested M_sqrt:", suggested, "\n")
  
  list(
    L_x = L_x,
    L_y = L_y,
    M_sqrt_needed = M_sqrt_needed,
    M_sqrt_suggested = suggested
  )
}


# =============================================================================
# USAGE EXAMPLE
# =============================================================================

# # 1. Run recommended workflow
result <- recommended_workflow(grid_coords_centered, grid_sf)
# 
# # 2. If happy with quantile trimming:
# new_params <- compute_new_L(result, c_factor = 1.5)
# 
# # 3. Update your grid for Stan model:
# grid_coords_trimmed <- grid_coords_centered[result$keep_idx, ]
# grid_sf_trimmed <- grid_sf[result$keep_idx, ]
# 
# # 4. Update Stan data:
# stan_data$grid_coords <- grid_coords_trimmed
# stan_data$L_x <- new_params$L_x
# stan_data$L_y <- new_params$L_y
# stan_data$M_sqrt <- new_params$M_sqrt_suggested
# stan_data$n_grids_total <- nrow(grid_coords_trimmed)
# 
# # 5. Need to also update all other data structures that reference grid cells!

cat("\n=== Script loaded successfully ===\n")
cat("Available functions:\n")
cat("  - trim_domain_quantiles()\n")
cat("  - trim_domain_density()\n")
cat("  - recommended_workflow()\n")
cat("  - compute_new_L()\n")