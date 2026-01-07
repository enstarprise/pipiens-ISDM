
# ============================================================================
# LOAD PACKAGES
# ============================================================================
# Set CRAN mirror
options(repos = c(CRAN = "https://cloud.r-project.org"))

required_packages <- c(
  # Data manipulation and I/O
  "tidyverse",    # Suite of data science packages
  "dplyr",        # Data manipulation
  "tidyr", # Data tidying
  "posterior",       
  "bayesplot",    # Date/time manipulation
  
  # Modeling
  "cmdstanr",       #cmdstanr
  
  # Utilities
  "optparse"       # Progress bars for apply functions
)

suppressPackageStartupMessages({
  library(cmdstanr)
  library(optparse)
  library(dplyr)
  library(bayesplot)
  library(posterior)
})

cat("\n==============================================\n")
cat(" STAN MODEL EXECUTION - SLURM JOB\n")
cat("==============================================\n\n")

# ----------------------------------------------------------------------------
# SLURM INFO
# ----------------------------------------------------------------------------
job_id  <- Sys.getenv("SLURM_JOB_ID")
n_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
n_tasks <- as.integer(Sys.getenv("SLURM_NTASKS"))          # Number of tasks (chains)
cpus_per_task <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))  # CPUs per task

cat("SLURM Job ID:", job_id, "\n")
cat("Available cores:", n_cores, "\n\n")

# ============================================================================
# MODEL SETUP AND EXECUTION
# ============================================================================
cat("Setting up model configuration...\n")

# Set defaults if not running under SLURM
if (is.na(job_id)) job_id <- "local_run"
if (is.na(n_tasks)) n_tasks <- 4
if (is.na(cpus_per_task)) cpus_per_task <- 4

cat("SLURM Job ID:", job_id, "\n")
cat("Number of tasks (chains):", n_tasks, "\n")
cat("CPUs per task (threads per chain):", cpus_per_task, "\n")
cat("Total cores requested:", n_tasks * cpus_per_task, "\n\n")

# Use SLURM settings for threading
n_chains <- n_tasks
threads_per_chain <- cpus_per_task

# Set environment variables for threading - use SLURM settings
Sys.setenv(STAN_NUM_THREADS = threads_per_chain)
Sys.setenv(OMP_NUM_THREADS = threads_per_chain)

# Also set for cmdstanr
options(mc.cores = threads_per_chain)

cat("Threading configuration:\n")
cat("  Chains:", n_chains, "\n")
cat("  Threads per chain:", threads_per_chain, "\n")
cat("  Total threads:", n_chains * threads_per_chain, "\n")
cat("  STAN_NUM_THREADS:", Sys.getenv("STAN_NUM_THREADS"), "\n")
cat("  OMP_NUM_THREADS:", Sys.getenv("OMP_NUM_THREADS"), "\n\n")
# ============================================================================
# COMPILE MODEL
# ============================================================================
cat("Compiling Stan model...\n")
compile_start <- Sys.time()

load("01_data/processedCovariates/1km/stan_data_init_1km_full.RData")

# Use the path to your model file
model_file <- "~/pipiens-ISDM/02_model/stan/marginalised_obs_model.stan"

# Check if model file exists
if (!file.exists(model_file)) {
  stop("Model file not found at: ", model_file)
}

# Compile with threading support
model <- cmdstan_model(
  model_file,
  cpp_options = list(
    stan_threads = TRUE,
    "PRECOMPILED_HEADERS" = FALSE  # May help with compilation issues
  ),
  stanc_options = list("O1"),  # Optimization level
  force_recompile = FALSE  # Set to TRUE if you need to recompile
)

compile_time <- difftime(Sys.time(), compile_start, units = "mins")
cat("Model compiled in", round(compile_time, 1), "minutes\n\n")

# ============================================================================
# RUN SAMPLING
# ============================================================================
cat("Starting MCMC sampling...\n")
cat("Configuration:\n")
cat("  Chains:", n_chains, "\n")
cat("  Parallel chains:", n_chains, "\n")
cat("  Threads per chain:", threads_per_chain, "\n")
cat("  Grainsize:", stan_data$grainsize, "\n")
cat("  Warmup iterations: 500\n")
cat("  Sampling iterations: 1000\n\n")

sampling_start <- Sys.time()

# Run the sampling
fit <- tryCatch({
  model$sample(
    data = stan_data,
    chains = n_chains,
    parallel_chains = n_chains,  # Run chains in parallel
    threads_per_chain = threads_per_chain,  # Threads for within-chain parallelization
    iter_warmup = 500,
    iter_sampling = 1000,
    refresh = 50,  # Print progress every 50 iterations
    adapt_delta = 0.95,
    max_treedepth = 12,
    init = init_list,  # Use your initial values list
    seed = 12345,
    show_messages = TRUE,
    validate_csv = FALSE  # Can speed up reading
  )
}, error = function(e) {
  cat("Error during sampling:", e$message, "\n")
  # Try with fewer parallel chains if failed
  if (n_chains > 1) {
    cat("Trying with fewer parallel chains...\n")
    model$sample(
      data = stan_data,
      chains = n_chains,
      parallel_chains = min(2, n_chains),  # Reduce parallel chains
      threads_per_chain = threads_per_chain,
      iter_warmup = 500,
      iter_sampling = 1000,
      refresh = 50,
      adapt_delta = 0.95,
      max_treedepth = 12,
      init = init_list,
      seed = 12345,
      show_messages = TRUE
    )
  } else {
    stop("Sampling failed completely")
  }
})

sampling_time <- difftime(Sys.time(), sampling_start, units = "mins")
cat("\nSampling completed in", round(sampling_time, 1), "minutes\n\n")

# ============================================================================
# SAVE RESULTS
# ============================================================================
cat("Saving results...\n")

# Create output directory with job ID
output_dir <- paste0("results_", job_id)
dir.create(output_dir, showWarnings = FALSE)

# Save the fit object
saveRDS(fit, file = file.path(output_dir, "stan_fit.rds"))

# Save draws in a more efficient format
draws <- fit$draws()
saveRDS(draws, file = file.path(output_dir, "draws.rds"))

# Save summary statistics
summary_df <- fit$summary()
write.csv(summary_df, file = file.path(output_dir, "summary.csv"), row.names = FALSE)

# ============================================================================
# DIAGNOSTICS
# ============================================================================
cat("Generating diagnostics...\n")

# Create diagnostics report
sink(file.path(output_dir, "diagnostics.txt"))
cat("Stan Model Diagnostics\n")
cat("=====================\n")
cat("Job ID:", job_id, "\n")
cat("Run completed:", Sys.time(), "\n")
cat("Data: 5k pixel subset\n")
cat("Chains:", n_chains, "\n")
cat("Threads per chain:", threads_per_chain, "\n")
cat("Total iterations: warmup = 500, sampling = 1000\n")
cat("Total runtime:", round(as.numeric(compile_time + sampling_time), 1), "minutes\n\n")

# Check diagnostics
diag <- fit$diagnostic_summary()
cat("Divergent transitions:", sum(unlist(diag$num_divergent)), "\n")
cat("Max tree depth hits:", sum(unlist(diag$num_max_treedepth)), "\n\n")

# Convergence diagnostics
cat("Convergence diagnostics:\n")
cat("Maximum R-hat:", round(max(summary_df$rhat, na.rm = TRUE), 3), "\n")
cat("Minimum bulk ESS:", round(min(summary_df$ess_bulk, na.rm = TRUE)), "\n")
cat("Minimum tail ESS:", round(min(summary_df$ess_tail, na.rm = TRUE)), "\n\n")

# Parameter summary
cat("Parameter statistics:\n")
cat("Total parameters:", nrow(summary_df), "\n")
cat("Parameters with R-hat > 1.05:", sum(summary_df$rhat > 1.05, na.rm = TRUE), "\n")
cat("Parameters with ESS < 100:", sum(summary_df$ess_bulk < 100, na.rm = TRUE), "\n")
sink()

# ============================================================================
# CREATE DIAGNOSTIC PLOTS
# ============================================================================
cat("Creating diagnostic plots...\n")

# Set up plotting
plot_file <- file.path(output_dir, "diagnostic_plots.png")
png(plot_file, width = 1200, height = 900)

# Layout for multiple plots
layout(matrix(c(1,2,3,4,5,6,7,8,9), nrow=3, byrow=TRUE))

# 1. R-hat histogram
hist(summary_df$rhat, breaks = 30, col = "lightblue",
     main = "R-hat Distribution", xlab = "R-hat")
abline(v = 1.05, col = "red", lwd = 2, lty = 2)
abline(v = 1.01, col = "orange", lwd = 1.5, lty = 3)

# 2. ESS histogram
hist(summary_df$ess_bulk, breaks = 30, col = "lightgreen",
     main = "Bulk ESS Distribution", xlab = "Effective Sample Size")
abline(v = 100, col = "red", lwd = 2, lty = 2)
abline(v = 400, col = "orange", lwd = 1.5, lty = 3)

# 3. Trace plot for key parameters (first few)
param_names <- summary_df$variable
# Look for key parameters
key_params <- param_names[grep("sigma|tau|alpha|beta", param_names)]
if (length(key_params) > 0) {
  key_params <- key_params[1:min(4, length(key_params))]
  
  for (i in seq_along(key_params)) {
    if (key_params[i] %in% param_names) {
      draws_mat <- as_draws_matrix(draws)
      if (key_params[i] %in% colnames(draws_mat)) {
        trace_data <- draws_mat[, key_params[i]]
        # Reshape for trace plot
        trace_matrix <- matrix(trace_data, nrow = 1000, ncol = n_chains)
        matplot(trace_matrix, type = "l", lty = 1,
                col = 1:n_chains, 
                main = paste("Trace:", key_params[i]),
                xlab = "Iteration", ylab = "Value",
                cex.main = 0.8)
        if (i == 1) {
          legend("topright", legend = paste("Chain", 1:n_chains),
                 col = 1:n_chains, lty = 1, cex = 0.6)
        }
      }
    }
  }
}

# 4. Timing per chain
if (!is.null(fit$time()$chains)) {
  chain_times <- sapply(fit$time()$chains, function(x) x$total)
  barplot(chain_times, col = rainbow(n_chains),
          main = "Runtime per Chain", ylab = "Seconds",
          names.arg = paste("Chain", 1:n_chains))
}

dev.off()

# ============================================================================
# SAVE TIMING INFORMATION
# ============================================================================
timing_info <- data.frame(
  job_id = job_id,
  n_cores_total = n_cores,
  n_chains = n_chains,
  threads_per_chain = threads_per_chain,
  grainsize = stan_data$grainsize,
  compile_time_min = as.numeric(compile_time),
  sampling_time_min = as.numeric(sampling_time),
  total_time_min = as.numeric(compile_time + sampling_time),
  iterations_warmup = 500,
  iterations_sampling = 1000,
  date = Sys.time()
)

write.csv(timing_info, file = file.path(output_dir, "timing_info.csv"), row.names = FALSE)

# ============================================================================
# FINAL SUMMARY
# ============================================================================
cat("\n" , strrep("=", 60), "\n")
cat("JOB COMPLETED SUCCESSFULLY\n")
cat(strrep("=", 60), "\n")
cat("Results saved to:", output_dir, "\n\n")
cat("Files created:\n")
cat("  stan_fit.rds          - Complete Stan fit object\n")
cat("  draws.rds             - Posterior draws\n")
cat("  summary.csv           - Parameter summary statistics\n")
cat("  diagnostics.txt       - Diagnostic report\n")
cat("  diagnostic_plots.png  - Diagnostic plots\n")
cat("  timing_info.csv       - Timing information\n\n")

cat("Key diagnostics:\n")
cat("  Max R-hat:", round(max(summary_df$rhat, na.rm = TRUE), 3), "\n")
cat("  Min ESS:", round(min(summary_df$ess_bulk, na.rm = TRUE)), "\n")
cat("  Divergent transitions:", sum(unlist(diag$num_divergent)), "\n\n")

cat("Total runtime:", round(as.numeric(compile_time + sampling_time), 1), "minutes\n")
cat("Completed at:", Sys.time(), "\n")

# ============================================================================
# OPTIONAL: CLEANUP
# ============================================================================
# Remove large objects if memory is tight
# rm(draws, draws_mat)
# gc()