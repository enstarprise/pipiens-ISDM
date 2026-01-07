library(cmdstanr)
library(posterior)
library(dplyr)
library(tidyr)



# Define grainsizes to test
grainsizes_to_test <- c(100, 200, 400, 800, 1000)
# Or for finer testing: 
# grainsizes_to_test <- seq(50, 500, by = 50)

# Create a directory for benchmarking results
if (!dir.exists("benchmark_results")) {
  dir.create("benchmark_results")
}

# Load/compile model once
load("02_model/stan/01_5kpixel_test/5k_pixelSubset.RData"); stan_data <- subset_stan_data
model <- cmdstan_model("02_model/stan/marginalised_obs_model.stan",
                       cpp_options = list(stan_threads = TRUE))

# Store benchmarking results
benchmark_results <- data.frame(
  grainsize = numeric(),
  total_time = numeric(),
  sampling_time = numeric(),
  warmup_time = numeric(),
  cpu_utilization = numeric(),
  ess_per_sec = numeric(),
  divergences = numeric(),
  max_rhat = numeric(),
  min_ess = numeric(),
  chains = numeric(),
  threads_per_chain = numeric(),
  date = as.character(character())
)

# Function to extract key diagnostics
# extract_diagnostics <- function(fit) {
#   summ <- fit$summary()
#   diag <- fit$diagnostic_summary()
#   
#   return(list(
#     max_rhat = max(summ$rhat, na.rm = TRUE),
#     min_ess_bulk = min(summ$ess_bulk, na.rm = TRUE),
#     min_ess_tail = min(summ$ess_tail, na.rm = TRUE),
#     divergences = if (!is.null(diag$num_divergent)) sum(diag$num_divergent) else NA
#   ))
# }

# Function to run benchmark for a single grainsize
# Revised benchmarking function with better error handling
run_benchmark <- function(grainsize, save_samples = FALSE) {
  cat("\n" , strrep("=", 60), "\n")
  cat("Testing grainsize:", grainsize, "\n")
  cat(strrep("=", 60), "\n")
  
  # Set grainsize in data
  stan_data$grainsize <- grainsize
  stan_data$N_multiplier <- 5
  
  # Start timer
  start_time <- Sys.time()
  
  tryCatch({
    # Run sampling - store result in a variable
    fit_result <- model$sample(
      data = stan_data,
      chains = 2,  # Use 2 chains for faster benchmarking
      parallel_chains = 2,
      threads_per_chain = 4,
      iter_warmup = 50,    # Reduced for benchmarking
      iter_sampling = 100,  # Reduced for benchmarking
      refresh = 25,
      adapt_delta = 0.95,
      max_treedepth = 12,
      init = init_fun,
      seed = 123
    )
    
    # Check if fit_result is valid
    if (!inherits(fit_result, "CmdStanMCMC")) {
      stop("Model fitting did not return a valid CmdStanMCMC object")
    }
    
    # Calculate times
    total_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    # Safely extract timing information
    timing_info <- tryCatch({
      fit_result$time()
    }, error = function(e) {
      cat("Warning: Could not extract timing info:", e$message, "\n")
      return(list(chains = list(
        list(warmup = NA, sampling = NA),
        list(warmup = NA, sampling = NA)
      )))
    })
    
    # Calculate warmup and sampling times
    if (all(!is.na(timing_info$chains))) {
      warmup_time <- sum(sapply(timing_info$chains, function(x) x$warmup), na.rm = TRUE)
      sampling_time <- sum(sapply(timing_info$chains, function(x) x$sampling), na.rm = TRUE)
    } else {
      warmup_time <- NA
      sampling_time <- NA
    }
    
    # Safely extract diagnostics
    diag_info <- tryCatch({
      extract_diagnostics(fit_result)
    }, error = function(e) {
      cat("Warning: Could not extract diagnostics:", e$message, "\n")
      return(list(
        max_rhat = NA,
        min_ess_bulk = NA,
        min_ess_tail = NA,
        divergences = NA
      ))
    })
    
    # Calculate ESS per second
    ess_per_sec <- NA
    if (!is.na(sampling_time) && sampling_time > 0) {
      tryCatch({
        summ <- fit_result$summary()
        if (nrow(summ) > 0 && "ess_bulk" %in% names(summ)) {
          # Use all parameters or key parameters
          valid_ess <- summ$ess_bulk[!is.na(summ$ess_bulk) & is.finite(summ$ess_bulk)]
          if (length(valid_ess) > 0) {
            ess_per_sec <- mean(valid_ess) / (sampling_time / 2)  # per chain (2 chains)
          }
        }
      }, error = function(e) {
        cat("Warning: Could not calculate ESS/sec:", e$message, "\n")
      })
    }
    
    # Calculate CPU utilization (estimate)
    cpu_utilization <- if (!is.na(total_time) && !is.na(sampling_time) && total_time > 0) {
      (sampling_time / total_time) * 100
    } else {
      NA
    }
    
    # Save samples if requested
    if (save_samples && inherits(fit_result, "CmdStanMCMC")) {
      tryCatch({
        samples_file <- paste0("benchmark_results/samples_grainsize_", grainsize, ".rds")
        saveRDS(fit_result$draws(), samples_file)
        cat("Samples saved to:", samples_file, "\n")
      }, error = function(e) {
        cat("Warning: Could not save samples:", e$message, "\n")
      })
    }
    
    # Return results
    result <- data.frame(
      grainsize = grainsize,
      total_time = total_time,
      sampling_time = sampling_time,
      warmup_time = warmup_time,
      cpu_utilization = cpu_utilization,
      ess_per_sec = ess_per_sec,
      divergences = ifelse(is.null(diag_info$divergences), NA, diag_info$divergences),
      max_rhat = diag_info$max_rhat,
      min_ess = diag_info$min_ess_bulk,
      chains = 2,
      threads_per_chain = 4,
      date = as.character(Sys.time()),
      status = "success"
    )
    
    cat("\nResults:\n")
    cat("Total time:", round(total_time, 1), "seconds\n")
    if (!is.na(sampling_time)) cat("Sampling time:", round(sampling_time, 1), "seconds\n")
    if (!is.na(warmup_time)) cat("Warmup time:", round(warmup_time, 1), "seconds\n")
    if (!is.na(cpu_utilization)) cat("CPU utilization:", round(cpu_utilization, 1), "%\n")
    if (!is.na(ess_per_sec)) cat("ESS/sec:", round(ess_per_sec, 2), "\n")
    if (!is.na(diag_info$divergences)) cat("Divergences:", diag_info$divergences, "\n")
    if (!is.na(diag_info$max_rhat)) cat("Max R-hat:", round(diag_info$max_rhat, 3), "\n")
    
    return(result)
    
  }, error = function(e) {
    cat("Error with grainsize", grainsize, ":", e$message, "\n")
    
    # Return error result
    error_result <- data.frame(
      grainsize = grainsize,
      total_time = as.numeric(difftime(Sys.time(), start_time, units = "secs")),
      sampling_time = NA,
      warmup_time = NA,
      cpu_utilization = NA,
      ess_per_sec = NA,
      divergences = NA,
      max_rhat = NA,
      min_ess = NA,
      chains = 2,
      threads_per_chain = 4,
      date = as.character(Sys.time()),
      status = paste("error:", e$message)
    )
    
    return(error_result)
  })
}

# Helper function to extract diagnostics with safety checks
extract_diagnostics <- function(fit) {
  tryCatch({
    # Check if fit has summary method
    if (!inherits(fit, "CmdStanMCMC")) {
      stop("Not a valid CmdStanMCMC object")
    }
    
    # Get summary
    summ <- fit$summary()
    
    # Get diagnostic summary
    diag_summ <- tryCatch({
      fit$diagnostic_summary()
    }, error = function(e) {
      # Return empty list if diagnostic summary fails
      list()
    })
    
    # Calculate diagnostics
    max_rhat <- if (!is.null(summ$rhat) && length(summ$rhat) > 0) {
      max(summ$rhat, na.rm = TRUE)
    } else {
      NA
    }
    
    min_ess_bulk <- if (!is.null(summ$ess_bulk) && length(summ$ess_bulk) > 0) {
      min(summ$ess_bulk, na.rm = TRUE)
    } else {
      NA
    }
    
    min_ess_tail <- if (!is.null(summ$ess_tail) && length(summ$ess_tail) > 0) {
      min(summ$ess_tail, na.rm = TRUE)
    } else {
      NA
    }
    
    # Get divergences
    divergences <- if (!is.null(diag_summ$num_divergent)) {
      if (is.list(diag_summ$num_divergent)) {
        sum(unlist(diag_summ$num_divergent), na.rm = TRUE)
      } else {
        sum(diag_summ$num_divergent, na.rm = TRUE)
      }
    } else {
      NA
    }
    
    return(list(
      max_rhat = max_rhat,
      min_ess_bulk = min_ess_bulk,
      min_ess_tail = min_ess_tail,
      divergences = divergences
    ))
    
  }, error = function(e) {
    cat("Error in extract_diagnostics:", e$message, "\n")
    return(list(
      max_rhat = NA,
      min_ess_bulk = NA,
      min_ess_tail = NA,
      divergences = NA
    ))
  })
}



# Minimal benchmarking - just timing
simple_benchmark <- function(grainsizes = c(100, 200, 400, 800, 1000)) {
  results <- data.frame()
  
  for (gs in grainsizes) {
    cat("\n=== Testing grainsize:", gs, "===\n")
    
    stan_data$grainsize <- gs
    
    # Time the compilation separately
    compile_start <- Sys.time()
    
    tryCatch({
      # Run with minimal output
      fit <- model$sample(
        data = stan_data,
        chains = 1,
        # parallel_chains = 2,
        threads_per_chain = 4,
        iter_warmup = 50,
        iter_sampling = 100,
        refresh = 10,  # = 0 fir No output for cleaner timing
        adapt_delta = 0.95,
        max_treedepth = 12,
        init = init_fun,
        seed = 123,
        show_messages = TRUE
      )
      
      total_time <- as.numeric(difftime(Sys.time(), compile_start, units = "secs"))
      
      # Simple timing extraction
      if (exists("fit") && !is.null(fit)) {
        tryCatch({
          timing <- fit$time()
          warmup_time <- if (!is.null(timing)) {
            sum(sapply(timing$chains, function(x) x$warmup))
          } else {
            NA
          }
          sampling_time <- if (!is.null(timing)) {
            sum(sapply(timing$chains, function(x) x$sampling))
          } else {
            NA
          }
        }, error = function(e) {
          warmup_time <- NA
          sampling_time <- NA
        })
      }
      
      results <- rbind(results, data.frame(
        grainsize = gs,
        total_time = total_time,
        warmup_time = warmup_time,
        sampling_time = sampling_time
      ))
      
      cat("Total time:", round(total_time, 1), "seconds\n")
      
    }, error = function(e) {
      cat("ERROR:", e$message, "\n")
      results <- rbind(results, data.frame(
        grainsize = gs,
        total_time = as.numeric(difftime(Sys.time(), compile_start, units = "secs")),
        warmup_time = NA,
        sampling_time = NA
      ))
    })
    
    # Clean up memory
    if (exists("fit")) rm(fit)
    gc()
    
    # Pause between runs
    Sys.sleep(5)
  }
  
  return(results)
}

# Run simple benchmark
simple_results <- simple_benchmark()
print(simple_results)

# Plot results
if (nrow(simple_results) > 1) {
  plot(simple_results$grainsize, simple_results$total_time,
       type = "b", col = "blue", lwd = 2,
       xlab = "Grainsize", ylab = "Total Time (seconds)",
       main = "Grainsize Benchmarking")
  
  # Add smooth curve
  try({
    loess_fit <- loess(total_time ~ grainsize, data = simple_results)
    grainsize_seq <- seq(min(simple_results$grainsize), max(simple_results$grainsize), length = 100)
    lines(grainsize_seq, predict(loess_fit, data.frame(grainsize = grainsize_seq)), 
          col = "red", lty = 2)
  })
}





# Run benchmarks for all grainsizes
for (gs in grainsizes_to_test) {
  result <- run_benchmark(gs, save_samples = FALSE)
  
  if (!is.null(result)) {
    benchmark_results <- rbind(benchmark_results, result)
    
    # Save intermediate results
    write.csv(benchmark_results, "benchmark_results/grainsize_results.csv", 
              row.names = FALSE)
    
    # Plot progress
    if (nrow(benchmark_results) > 1) {
      png("benchmark_results/progress_plot.png", width = 1200, height = 800)
      par(mfrow = c(2, 2))
      
      # Plot 1: Time vs grainsize
      plot(benchmark_results$grainsize, benchmark_results$total_time,
           type = "b", col = "blue", lwd = 2,
           xlab = "Grainsize", ylab = "Total Time (s)",
           main = "Total Runtime vs Grainsize")
      
      # Plot 2: ESS/sec vs grainsize
      plot(benchmark_results$grainsize, benchmark_results$ess_per_sec,
           type = "b", col = "green", lwd = 2,
           xlab = "Grainsize", ylab = "ESS per second",
           main = "Efficiency vs Grainsize")
      
      # Plot 3: CPU utilization vs grainsize
      plot(benchmark_results$grainsize, benchmark_results$cpu_utilization,
           type = "b", col = "red", lwd = 2,
           xlab = "Grainsize", ylab = "CPU Utilization (%)",
           main = "CPU Utilization vs Grainsize",
           ylim = c(0, 100))
      
      # Plot 4: Sampling time vs grainsize
      plot(benchmark_results$grainsize, benchmark_results$sampling_time,
           type = "b", col = "purple", lwd = 2,
           xlab = "Grainsize", ylab = "Sampling Time (s)",
           main = "Sampling Time vs Grainsize")
      
      dev.off()
    }
  }
  
  # Optional: Add a cooldown period between runs
  Sys.sleep(10)
}

# Print final summary
cat("\n" , strrep("=", 60), "\n")
cat("BENCHMARKING COMPLETE\n")
cat(strrep("=", 60), "\n\n")

print(benchmark_results)

# Find optimal grainsize based on ESS/sec
if (nrow(benchmark_results) > 0) {
  optimal_idx <- which.max(benchmark_results$ess_per_sec)
  optimal_gs <- benchmark_results$grainsize[optimal_idx]
  
  cat("\nOptimal grainsize:", optimal_gs, "\n")
  cat("Maximum ESS/sec:", round(benchmark_results$ess_per_sec[optimal_idx], 2), "\n")
  cat("Total time at optimal:", round(benchmark_results$total_time[optimal_idx], 1), "seconds\n")
  
  # Generate final report
  sink("benchmark_results/final_report.txt")
  cat("Grainsize Benchmarking Report\n")
  cat("Generated:", Sys.time(), "\n\n")
  cat("Model: marginalised_obs_model.stan\n")
  cat("Data dimensions: 85,603 pixels × 178 dates\n")
  cat("Total observations:", 85603 * 178, "\n\n")
  cat("Optimal grainsize:", optimal_gs, "\n\n")
  cat("All results:\n")
  print(benchmark_results)
  sink()
}

# Alternative: Quick test with fewer iterations (for faster screening)
quick_benchmark <- function(grainsizes = c(50, 100, 200, 400)) {
  quick_results <- data.frame()
  
  for (gs in grainsizes) {
    cat("\nQuick test - grainsize:", gs, "\n")
    
    stan_data$grainsize <- gs
    stan_data$N_multiplier <- 5
    
    start <- Sys.time()
    
    fit_quick <- model$sample(
      data = stan_data,
      chains = 2,
      parallel_chains = 2,
      threads_per_chain = 4,
      iter_warmup = 50,
      iter_sampling = 100,
      refresh = 20,
      adapt_delta = 0.95,
      max_treedepth = 10,
      init = init_fun,
      seed = 123
    )
    
    total_time <- as.numeric(difftime(Sys.time(), start, units = "secs"))
    
    quick_results <- rbind(quick_results, data.frame(
      grainsize = gs,
      total_time = total_time
    ))
  }
  
  return(quick_results)
}

# Uncomment to run quick benchmark first:
# quick_results <- quick_benchmark()
# print(quick_results)