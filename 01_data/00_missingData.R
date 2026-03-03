###############################################################################
# MISSING VALUES
# see 2nd edition of statistical rethinking
# https://mc-stan.org/docs/stan-users-guide/missing-data.html
# https://cran.r-project.org/web/packages/finalfit/vignettes/missing.html
# https://stefvanbuuren.name/fimd/missing-data-pattern.html

# the covariates that have missing values are temperature, rainfall, livestock,
# elevation

library(finalfit)
###############################################################################
load("01_data/processedCovariates/1km/loaded_data_1km.RData")
load("01_data/processedCovariates/1km/env_land_data_1km.RData")
load("01_data/processedCovariates/1km/env_climate_data_1km.RData")
load("01_data/processedCovariates/1km/env_temp_data_1km.RData")
load("01_data/processedCovariates/1km/env_rain_data_1km_part1.RData")
load("01_data/processedCovariates/1km/env_rain_data_1km_part2.RData")


# ============================================================================
# 1. Convert to long format dataframe for finalfit
# ============================================================================
z_land <- as.data.frame(z_land_data) %>%
  dplyr::arrange(grid_id) %>% # grid_id_area2km
  dplyr::select(starts_with(c("grid_id", "z_"))) %>%
  dplyr::select(-c(z_saltwater_km2, z_urban_km2, z_suburban_km2, z_poi_count, z_reports)) %>% # z_poi_count_grouped
  as.matrix()

df_land <- as.data.frame(z_land)

# For temporal variables (z_temp)
# Assuming z_temp dimensions are: [grid, time, variable]
z_climate_data <- z_climate_data[,1:4]
n_grids <- dim(z_climate_data)[1]
n_times <- dim(z_climate_data)[2]
n_vars <- 2  # Two climate variables in columns 3:4

colnames(z_climate_data) <- c("grid_id", "date", "temperature", "precipitation")

df_climate_long <- as.data.frame(z_climate_data) %>%
  pivot_longer(
    cols = c(temperature, precipitation),
    names_to = "variable",
    values_to = "value"
  )


df_full <- z_climate_data %>%
  left_join(df_land, by = "grid_id")



# ============================================================================
# 2. Missing data patterns (Statistical Rethinking approach)
# ============================================================================

# Key question from McElreath: Is missingness related to the value itself (MNAR)
# or to other observed variables (MAR)?

# Create missing indicators for each variable
vars_to_check <- setdiff(names(df_full), c("date"))
for(var in vars_to_check) {
  df_full[[paste0(var, "_missing")]] <- is.na(df_full[[var]])
}

# ============================================================================
# 3. Visualize missing patterns with finalfit
# ============================================================================

# Overall missing pattern
missing_pattern <- df_full %>%
  select(all_of(vars_to_check)) %>%
  ff_glimpse()

# Missing pattern plot
missing_plot <- df_full %>%
  select(all_of(vars_to_check)) %>%
  missing_pattern()

print(missing_plot)

# ============================================================================
# 4. Test for MAR vs MCAR (Missing Completely at Random)
# ============================================================================
# Compare observed vs missing groups for each variable
# If other variables differ between missing/non-missing groups, suggests MAR

missing_compare_list <- list()
for(var in vars_to_check) {
  missing_indicator <- paste0(var, "_missing")
  
  # Skip if no missingness
  if(sum(df_full[[missing_indicator]]) == 0 || 
     sum(df_full[[missing_indicator]]) == nrow(df_full)) {
    missing_compare_list[[var]] <- "No missingness or completely missing"
    next
  }
  
  # Compare all other variables by missing status
  explanatory <- setdiff(vars_to_check, var)
  
  # For continuous variables, use t-tests/wilcoxon
  comparison_df <- data.frame(
    variable = character(),
    missing_mean = numeric(),
    observed_mean = numeric(),
    p_value = numeric(),
    test = character(),
    stringsAsFactors = FALSE
  )
  
  for(exp_var in explanatory) {
    # Skip if explanatory variable has too much missingness
    if(sum(is.na(df_full[[exp_var]])) > 0.9 * nrow(df_full)) next
    
    # Remove rows where explanatory variable is also missing
    temp_df <- df_full[!is.na(df_full[[exp_var]]), ]
    
    if(nrow(temp_df) < 10) next
    
    # Check if variable is numeric
    if(is.numeric(temp_df[[exp_var]])) {
      missing_vals <- temp_df[[exp_var]][temp_df[[missing_indicator]]]
      observed_vals <- temp_df[[exp_var]][!temp_df[[missing_indicator]]]
      
      # Use Wilcoxon test (non-parametric, safer)
      tryCatch({
        test_result <- wilcox.test(missing_vals, observed_vals)
        
        comparison_df <- rbind(comparison_df, data.frame(
          variable = exp_var,
          missing_mean = mean(missing_vals, na.rm = TRUE),
          observed_mean = mean(observed_vals, na.rm = TRUE),
          p_value = test_result$p.value,
          test = "Wilcoxon",
          stringsAsFactors = FALSE
        ))
      }, error = function(e) NULL)
    }
  }
  
  # Sort by p-value
  if(nrow(comparison_df) > 0) {
    comparison_df <- comparison_df %>% arrange(p_value)
    missing_compare_list[[var]] <- comparison_df
  } else {
    missing_compare_list[[var]] <- "Could not perform comparisons"
  }
}

# Print summary of MAR tests
cat("\n=== TESTING FOR MAR (Missing at Random) ===\n")
for(var in names(missing_compare_list)) {
  cat("\n", var, ":\n")
  result <- missing_compare_list[[var]]
  if(is.character(result)) {
    cat("  ", result, "\n")
  } else if(is.data.frame(result)) {
    sig_results <- result[result$p_value < 0.05, ]
    if(nrow(sig_results) > 0) {
      cat("  *** Significant predictors of missingness (p < 0.05):\n")
      print(sig_results[, c("variable", "missing_mean", "observed_mean", "p_value")])
      cat("  -> Suggests MAR (missingness related to other variables)\n")
    } else {
      cat("  No significant predictors found -> Suggests MCAR\n")
    }
  }
}

save(missing_compare_list, comparison_df, comparison, sig_results, missing_plot, missing_pattern,
     file = "01_data/descriptivePlots/missingness_all_results.RData")
# save(missing_compare_list, file ="01_data/descriptivePlots/missing_compare_list.RData")
# save(comparison_df, file ="01_data/descriptivePlots/missing_comparison_df.RData")
# save(comparison, file ="01_data/descriptivePlots/missing_comparison.RData")
# save(sig_results, file ="01_data/descriptivePlots/missing_sig_results.RData")
# save(missing_pattern, file ="01_data/descriptivePlots/missing_pattern.RData")

# ============================================================================
# 5. Spatial pattern of missingness
# ============================================================================
library(sf)
library(dplyr)

grid_sf <- st_read("01_data/grids/grid_clipped_1km.gpkg")

df_full <- df_full %>%
  left_join(grid_coords, by = "grid_id")


# Aggregate missingness by grid
missing_by_grid <- df_full %>%
  group_by(grid_id) %>%
  summarise(across(all_of(vars_to_check), 
                   ~mean(is.na(.)), 
                   .names = "prop_missing_{.col}"))

# Visualize if you have spatial coordinates
# If z_land contains lat/lon:
if("latitude" %in% names(df_land) && "longitude" %in% names(df_land)) {
  missing_by_grid <- missing_by_grid %>%
    left_join(df_land %>% select(grid_id, latitude, longitude), 
              by = "grid_id")
  
  # Plot spatial pattern of missingness for a specific variable
  # Example for first variable:
  first_var <- vars_to_check[1]
  ggplot(missing_by_grid, 
         aes(x = longitude, y = latitude, 
             fill = get(paste0("prop_missing_", first_var)))) +
    geom_tile() +
    scale_fill_viridis_c(name = "Proportion\nMissing") +
    labs(title = paste("Spatial Pattern of Missingness:", first_var)) +
    theme_minimal()
}

# ============================================================================
# 6. Temporal pattern of missingness
# ============================================================================

missing_by_time <- df_full %>%
  group_by(date) %>%
  summarise(across(all_of(vars_to_check), 
                   ~mean(is.na(.)), 
                   .names = "prop_missing_{.col}"))

# Plot temporal pattern
missing_by_time %>%
  pivot_longer(cols = starts_with("prop_missing_"),
               names_to = "variable",
               values_to = "prop_missing") %>%
  mutate(variable = gsub("prop_missing_", "", variable)) %>%
  ggplot(aes(x = date, y = prop_missing, color = variable, group = variable)) +
  geom_line() +
  geom_point() +
  labs(title = "Temporal Pattern of Missingness",
       y = "Proportion Missing",
       x = "Time") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ============================================================================
# 7. Statistical Rethinking decision framework
# ============================================================================

# Based on McElreath's guidance:
# 
# Option 1: MCAR (Missing Completely at Random)
#   - Missingness unrelated to any variables
#   - Decision: Complete case analysis OK, or impute with grand mean
#
# Option 2: MAR (Missing at Random)
#   - Missingness related to OTHER observed variables
#   - Decision: Multiple imputation or model missingness explicitly
#   - In Stan/brms: impute during modeling with priors
#
# Option 3: MNAR (Missing Not at Random)
#   - Missingness related to the unobserved value itself
#   - Decision: Need domain knowledge; may need selection model
#   - Example: High values systematically missing (ceiling effects)

# To decide, examine:
cat("\n=== DECISION GUIDE ===\n")
cat("1. Check missing_compare_list: Do other variables predict missingness?\n")
cat("2. Check spatial/temporal patterns: Is missingness clustered?\n")
cat("3. Domain knowledge: Why might data be missing?\n\n")

cat("If MAR: Consider these approaches for your Bayesian model:\n")
cat("  - Multiple imputation pre-model (mice package)\n")
cat("  - Joint modeling in Stan (impute within model)\n")
cat("  - Include missing indicators as predictors\n\n")

cat("If MNAR: May need to:\n")
cat("  - Model the missing data mechanism explicitly\n")
cat("  - Use sensitivity analysis\n")
cat("  - Consult Statistical Rethinking Ch 15 for selection models\n")

# ============================================================================
# 8. Export results for further analysis
# ============================================================================

results <- list(
  missing_pattern = missing_pattern,
  missing_compare = missing_compare_list,
  # missing_by_grid = missing_by_grid,
  # missing_by_time = missing_by_time,
  df_with_indicators = df_full
)

results
