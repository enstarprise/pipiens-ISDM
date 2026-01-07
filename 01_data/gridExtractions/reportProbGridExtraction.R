library(raster)
library(terra) 
library(sf)
library(tidyr)

grid_clipped <- st_read("grids/grid_clipped_1km.gpkg")
grid_clipped <- st_transform(grid_clipped, st_crs(scotland))
st_crs(grid_clipped)

citizenScience_all <- read_csv("~/Library/CloudStorage/OneDrive-UniversityofGlasgow/PhD/modelling/covariateExtraction/covariate extraction/citizenScience_data_humanpop.csv")
citizenScience_all$Date_found <- as.Date(citizenScience_all$Date_found, 
                                         format = "%d/%m/%Y")
citizenScience_all$Year <- format(citizenScience_all$Date_found, "%Y")
citizenScience_all$Month <- format(citizenScience_all$Date_found, "%m")
citizenScience_all$presence <- 1 # one obs = 1 report

# some lat/lon missing
citizenScience_all <- citizenScience_all %>%
  drop_na(c("Longitude", "Latitude")) 


cs_all_sf <- st_as_sf(citizenScience_all, coords = c("Longitude", "Latitude"), crs = 4326)
class(cs_all_sf)
st_crs(cs_all_sf)

grid_sf <- st_as_sf(grid_clipped)
grid_sf <- st_transform(grid_sf, st_crs(cs_all_sf))
st_crs(grid_sf)

# count reports per grid cell (handling empty cells)
joined <- st_join(grid_sf, cs_all_sf, join = st_within)
grid_counts <- joined %>%
  group_by(grid_id) %>%    
  summarise(reports = sum(presence, na.rm = TRUE))  # sums 1s, 0 if empty
grid_counts <- as.data.frame(grid_counts)

grid_1km_reports <- grid_clipped
grid_1km_reports <- grid_1km_reports %>%
  left_join(
    grid_counts,
    by = "grid_id"
  ) 
# %>%
  # mutate(
  #   reports = coalesce(reports, 0) 
  # )

write_csv(grid_1km_reports, "csvs/CSreports_1km.csv")
