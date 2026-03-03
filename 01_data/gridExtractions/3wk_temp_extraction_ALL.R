# =============================================================================
# TEMPERATURE EXTRACTION — tasmin, tasmax, tasmean
# 21-day rolling mean for each variable across Scotland 1km grid
# see 3WK_... R scripts for the full workflow incl. examples and 
# sanity/data checks
# =============================================================================

# ---- CONFIG ------------------------------------------------------------------

nc_dir         <- "01_data/gridExtractions/covariates/daily/"
output_dir     <- "01_data/csvs/"
grid_path      <- "01_data/grids/grid_clipped_1km.gpkg"
scotland_path  <- "01_data/gridExtractions/scotlandBoundary/scotland_boundary.shp"
survey_path    <- "01_data/csvs/survey_df.csv"
cs_path        <- "01_data/csvs/cs_df.csv"

window_days    <- 21          #  window size
n_cores        <- detectCores() - 4

extract_tasmin <- TRUE        # extract minimum temperature?
extract_tasmax <- TRUE        # extract maximum temperature?
calc_tasmean   <- TRUE        # calculate & save combined mean? (requires both above)

# ------------------------------------------------------------------------------

library(sf)
library(ncdf4)
library(terra)
library(future.apply)
library(parallel)
library(pbapply)
library(abind)
library(lubridate)
library(stringr)
library(dplyr)
library(readr)


# ---- FUNCTIONS ---------------------------------------------------------------

load_nc_temp_data <- function(files, var_name) {
  lapply(files, function(file) {
    nc <- nc_open(file)
    temperature <- ncvar_get(nc, var_name)
    time        <- ncvar_get(nc, "time")
    units       <- ncatt_get(nc, "time", "units")$value
    origin      <- as.POSIXct(sub("hours since ", "", units), tz = "UTC")
    dates       <- as.Date(origin + hours(time))
    nc_close(nc)
    list(file = file, temperature = temperature, dates = dates)
  })
}

find_nearest_grid_bng <- function(x, y, x_coords, y_coords) {
  c(which.min(abs(x_coords - x)),
    which.min(abs(y_coords - y)))
}

filter_nc_files <- function(nc_dir, var_name, start_date, end_date) {
  pattern    <- paste0(var_name, "_hadukgrid_uk_1km_day_.*\\.nc$")
  nc_files   <- list.files(nc_dir, pattern = pattern, full.names = TRUE)
  file_dates <- str_extract(nc_files, "\\d{8}-\\d{8}")
  
  keep <- sapply(seq_along(nc_files), function(i) {
    if (is.na(file_dates[i])) return(FALSE)
    parts      <- strsplit(file_dates[i], "-")[[1]]
    file_start <- as.Date(parts[1], format = "%Y%m%d")
    file_end   <- as.Date(parts[2], format = "%Y%m%d")
    (file_start <= end_date) && (file_end >= start_date)
  })
  
  nc_files[keep]
}

get_unique_extraction_dates <- function(..., date_columns) {
  dfs  <- list(...)
  all_dates <- unlist(mapply(function(df, col) df[[col]], dfs, date_columns,
                             SIMPLIFY = FALSE))
  sort(unique(as.Date(all_dates)))
}

extract_temp_grid <- function(grid_df, dates_to_extract, nc_dir, var_name,
                              window_days, col_name, n_cores) {
  lookback <- window_days - 1
  
  # Spatial reference from first file
  nc_ref   <- nc_open(list.files(nc_dir,
                                 pattern = paste0(var_name, ".*\\.nc$"),
                                 full.names = TRUE)[1])
  x_coords <- ncvar_get(nc_ref, "projection_x_coordinate")
  y_coords <- ncvar_get(nc_ref, "projection_y_coordinate")
  nc_close(nc_ref)
  
  # Pre-compute nearest grid indices for all cells
  grid_idx <- t(apply(grid_df[, c("x_bng", "y_bng")], 1, function(coord) {
    find_nearest_grid_bng(coord[1], coord[2], x_coords, y_coords)
  }))
  colnames(grid_idx) <- c("grid_x", "grid_y")
  grid_df <- cbind(grid_df, grid_idx)
  
  cat(sprintf("[%s] Processing %d dates across %d grid cells with %d cores\n",
              var_name, length(dates_to_extract), nrow(grid_df), n_cores))
  
  cl <- makeCluster(n_cores)
  clusterExport(cl, c("grid_df", "nc_dir", "var_name", "lookback",
                      "filter_nc_files", "load_nc_temp_data",
                      "find_nearest_grid_bng", "col_name"),
                envir = environment())
  clusterEvalQ(cl, {
    library(ncdf4); library(abind); library(stringr); library(lubridate)
  })
  
  result <- pblapply(dates_to_extract, function(target_date) {
    start_date <- target_date - lookback
    nc_files   <- filter_nc_files(nc_dir, var_name, start_date, target_date)
    if (length(nc_files) == 0) return(NULL)
    
    temp_data  <- load_nc_temp_data(nc_files, var_name)
    all_dates  <- unlist(lapply(temp_data, `[[`, "dates"))
    temp_stack <- abind::abind(lapply(temp_data, `[[`, "temperature"), along = 3)
    date_mask  <- (all_dates >= start_date) & (all_dates <= target_date)
    if (!any(date_mask)) return(NULL)
    
    valid_mask  <- !is.na(grid_df$grid_x) & !is.na(grid_df$grid_y)
    temp_values <- rep(NA_real_, nrow(grid_df))
    
    if (any(valid_mask)) {
      indices <- cbind(grid_df$grid_x[valid_mask], grid_df$grid_y[valid_mask])
      temp_values[valid_mask] <- apply(indices, 1, function(idx) {
        mean(temp_stack[idx[1], idx[2], date_mask], na.rm = TRUE)
      })
    }
    
    data.frame(
      grid_id    = grid_df$grid_id,
      date       = target_date,
      x_bng      = grid_df$x_bng,
      y_bng      = grid_df$y_bng,
      value      = round(temp_values, 1),
      start_date = start_date,
      end_date   = target_date,
      stringsAsFactors = FALSE
    )
  }, cl = cl)
  
  stopCluster(cl)
  
  out <- do.call(rbind, result[!sapply(result, is.null)])
  names(out)[names(out) == "value"] <- col_name
  out
}


# ---- LOAD DATA ---------------------------------------------------------------

scotland   <- st_read(scotland_path) |> st_transform(crs = 27700)
grid_clipped <- st_read(grid_path)

coords   <- st_coordinates(st_centroid(grid_clipped))
grid_df  <- data.frame(
  grid_id = grid_clipped$grid_id,
  x_bng   = coords[, 1],
  y_bng   = coords[, 2]
)

survey_df <- read_csv(survey_path)
survey_df$Setup_date <- as.Date(survey_df$Setup_date, format = "%d/%m/%Y")

cs_df <- read_csv(cs_path) |>
  filter(Verified_mosquito == "Yes", Species == "pipiens", Stage == "Adult")
cs_df$Date_found <- as.Date(cs_df$Date_found, format = "%d/%m/%Y")

dates_to_extract <- get_unique_extraction_dates(
  survey_df, cs_df,
  date_columns = c("Setup_date", "Date_found")
)

cat(sprintf("Dates to extract: %d (%s to %s)\n",
            length(dates_to_extract),
            min(dates_to_extract), max(dates_to_extract)))


# ---- EXTRACTION --------------------------------------------------------------

col_suffix <- paste0(window_days, "d_celsius")

if (extract_tasmin) {
  tasmin_grid <- extract_temp_grid(
    grid_df         = grid_df,
    dates_to_extract = dates_to_extract,
    nc_dir          = nc_dir,
    var_name        = "tasmin",
    window_days     = window_days,
    col_name        = paste0("mean_min_temp_", col_suffix),
    n_cores         = n_cores
  )
  write.csv(tasmin_grid,
            file      = file.path(output_dir, paste0("tasmin_grid_", col_suffix, ".csv")),
            row.names = FALSE)
  cat(sprintf("[tasmin] Saved %d rows\n", nrow(tasmin_grid)))
}

if (extract_tasmax) {
  tasmax_grid <- extract_temp_grid(
    grid_df          = grid_df,
    dates_to_extract = dates_to_extract,
    nc_dir           = nc_dir,
    var_name         = "tasmax",
    window_days      = window_days,
    col_name         = paste0("mean_max_temp_", col_suffix),
    n_cores          = n_cores
  )
  write.csv(tasmax_grid,
            file      = file.path(output_dir, paste0("tasmax_grid_", col_suffix, ".csv")),
            row.names = FALSE)
  cat(sprintf("[tasmax] Saved %d rows\n", nrow(tasmax_grid)))
}


# ---- TASMEAN ------------------------------------------------------

if (calc_tasmean) {
  if (!extract_tasmin || !extract_tasmax) {
    # Try loading from disk if not extracted in this run
    tasmin_path <- file.path(output_dir, paste0("tasmin_grid_", col_suffix, ".csv"))
    tasmax_path <- file.path(output_dir, paste0("tasmax_grid_", col_suffix, ".csv"))
    
    if (!exists("tasmin_grid")) {
      stopifnot("tasmin CSV not found — run with extract_tasmin = TRUE first" = file.exists(tasmin_path))
      tasmin_grid <- read_csv(tasmin_path)
    }
    if (!exists("tasmax_grid")) {
      stopifnot("tasmax CSV not found — run with extract_tasmax = TRUE first" = file.exists(tasmax_path))
      tasmax_grid <- read_csv(tasmax_path)
    }
  }
  
  min_col  <- paste0("mean_min_temp_", col_suffix)
  max_col  <- paste0("mean_max_temp_", col_suffix)
  mean_col <- paste0("mean_temp_",     col_suffix)
  
  tasmean_grid <- tasmax_grid[, c("grid_id", "date", max_col)] |>
    inner_join(tasmin_grid[, c("grid_id", "date", min_col)],
               by = c("grid_id", "date")) |>
    mutate(!!mean_col := round((.data[[max_col]] + .data[[min_col]]) / 2, 1))
  
  write.csv(tasmean_grid,
            file      = file.path(output_dir, paste0("tasmean_grid_", col_suffix, ".csv")),
            row.names = FALSE)
  cat(sprintf("[tasmean] Saved %d rows\n", nrow(tasmean_grid)))
}


