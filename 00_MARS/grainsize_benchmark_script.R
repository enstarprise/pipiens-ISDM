#!/usr/bin/env Rscript

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
cat("SLURM Job ID:", job_id, "\n")
cat("Available cores:", n_cores, "\n\n")

# ----------------------------------------------------------------------------
# CLI OPTIONS
# ----------------------------------------------------------------------------
option_list <- list(
  make_option("--mode", default="run"),
  make_option("--grainsize", type="integer", default=500),
  make_option("--threads-per-chain", type="integer", default=4)
)

opt     <- parse_args(OptionParser(option_list = option_list))
mode    <- opt$mode
gs      <- opt$grainsize
threads <- opt$`threads-per-chain`

cat("Mode:", mode, "\n")
cat("Grainsize:", gs, "\n")
cat("Threads per chain:", threads, "\n\n")

# ----------------------------------------------------------------------------
# LOAD DATA (already subset to 5k pixels)
# ----------------------------------------------------------------------------
load("~/pipiens-ISDM/stan/5kpixel_stan_data.RData")

stan_data$grainsize <- gs

cat("Stan data loaded successfully\n")
cat("  n_grids_total:", stan_data$n_grids_total, "\n")
cat("  n_dates:", stan_data$n_dates, "\n")
cat("  n_obs_y:", stan_data$n_obs_y, "\n\n")

# Generate initial values
n_chains <- min(4, n_cores)
init_list <- replicate(n_chains, init_fun(), simplify = FALSE)

# ----------------------------------------------------------------------------
# MODEL COMPILE
# ----------------------------------------------------------------------------

stan_file <- "~/pipiens-ISDM/stan/stan_model.stan"

stopifnot("Stan file not found!" = file.exists(stan_file))

model <- cmdstan_model(stan_file, cpp_options = list(stan_threads = TRUE))
cat("Model compiled ️\n\n")


# ============================================================================
# BENCHMARK MODE
# ============================================================================
if (mode == "benchmark") {
  cat("\n Benchmarking grainsize:", gs, "...\n")
  
  timing <- system.time(
    model$sample(
      data = stan_data,
      chains = 1,
      parallel_chains = 1,
      threads_per_chain = threads,
      iter_warmup = 50,
      iter_sampling = 100,
      refresh = 10
    )
  )[["elapsed"]]
  
  cat("Result:", timing, "seconds\n")
  
  # append to single CSV
  dir.create("~/pipiens-ISDM/stan/benchmark", FALSE, FALSE)
  out_file <- "~/pipiens-ISDM/stan/benchmark/grainsize_benchmark.csv"
  
  if (!file.exists(out_file)) write("grainsize,seconds", out_file)
  write(sprintf("%s,%s", gs, timing), out_file, append = TRUE)
  
  cat("Saved:", out_file, "\n\n")
  quit(save = "no")
}

# ============================================================================
# FULL MODEL RUN
# ============================================================================
cat("Starting MCMC...\n\n")

fit <- model$sample(
  data = stan_data,
  init = init_list,
  chains = n_chains,
  parallel_chains = n_chains,
  threads_per_chain = threads,
  iter_warmup = 200,
  iter_sampling = 500,
  adapt_delta = 0.9,
  max_treedepth = 12,
  refresh = 50,
  seed = 123
)

fit$save_object("fit_lastrun.rds")
cat("\n Model finished. Results saved to fit_lastrun.rds ✨\n")
