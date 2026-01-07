
library(sf)
library(ggplot2)
library(dplyr)
library(mapview)
library(terra)
library(ggspatial)


grid_clipped <- st_read("grids/grid_clipped_2point5km.gpkg")
scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp")
crs(scotland)
scotland <- st_transform(scotland,crs = 27700)

###### --- grid map, 2.5km pixels ---  ######
bbox <- st_bbox(scotland) # or expected region
grid_clipped <- st_crop(grid_clipped, bbox)


ggplot() +
  geom_sf(data = grid_clipped, linewidth = 0.1) +
  theme_bw() +
  coord_sf(expand = FALSE)


bbox <- st_bbox(grid_clipped)
# choose how much to crop (e.g. 20% of the width from the left)
xrange <- bbox["xmax"] - bbox["xmin"]
yrange <- bbox["ymax"] - bbox["ymin"]
crop_fraction_left <- 0.35  # increase to crop more
crop_fraction_top <- 0.09

grid_plot <- ggplot() +
  geom_sf(data = grid_clipped, linewidth = 0.1) +
  geom_sf(data = scotland, fill = NA, color = "black", linewidth = 0.3) +
  theme_bw() +
  coord_sf(
    xlim = c(bbox["xmin"] + xrange * crop_fraction_left, bbox["xmax"]),
    ylim = c(bbox["ymin"], bbox["ymax"] - yrange * crop_fraction_top)
  ) +
  # add a scale bar (units depend on CRS; CRS=27700 uses meters)
  annotation_scale(
    location = "bl",       # bottom-left corner
    width_hint = 0.3,      # controls size
    text_cex = 0.8
  ) +
  # add a north arrow
  annotation_north_arrow(
    location = "br",       # bottom-right corner
    which_north = "true",
    style = north_arrow_fancy_orienteering,
    height = unit(1, "cm"),
    width = unit(1, "cm")
  )

grid_plot <- grid_plot +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background  = element_rect(fill = "transparent", colour = NA)
  )

ggsave(filename = "grid_plot.png", 
       plot = grid_plot, width = 6, height = 6, 
       units = "in", dpi = 300, bg = "transparent")


##### --- zoom in on some cells --- #####
bbox <- st_bbox(scotland)
grid_clipped <- st_crop(grid_clipped, bbox)

# Find the center of the grid
grid_center <- st_centroid(st_union(grid_clipped))

# Create a buffer around the center (adjust distance as needed)
zoom_area <- st_buffer(grid_center, dist = 25000)  # 25km buffer

# Crop grid to the zoom area
grid_zoomed <- st_intersection(grid_clipped, zoom_area)

# Plot the zoomed area
ggplot() +
  geom_sf(data = grid_zoomed, linewidth = 0.1) +
  geom_sf(data = scotland, fill = NA, color = "black", linewidth = 0.3) +
  theme_bw() +
  coord_sf(expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.8) +
  annotation_north_arrow(location = "br", which_north = "true",
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"), width = unit(1, "cm"))


create_zoomed_plot <- function(zoom_factor = 0.5) {
  full_bbox <- st_bbox(grid_clipped)
  xrange <- full_bbox["xmax"] - full_bbox["xmin"]
  yrange <- full_bbox["ymax"] - full_bbox["ymin"]
  
  x_crop <- (1 - zoom_factor) / 2
  y_crop <- (1 - zoom_factor) / 2
  
  new_xlim <- c(full_bbox["xmin"] + xrange * x_crop, 
                full_bbox["xmax"] - xrange * x_crop)
  new_ylim <- c(full_bbox["ymin"] + yrange * y_crop, 
                full_bbox["ymax"] - yrange * y_crop)
  
  ggplot() +
    geom_sf(data = grid_clipped, linewidth = 0.1) +
    geom_sf(data = scotland, fill = NA, color = "black", linewidth = 0.3) +
    theme_bw() +
    coord_sf(xlim = new_xlim, ylim = new_ylim, expand = FALSE) +
    annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.8) +
    annotation_north_arrow(location = "br", which_north = "true",
                           style = north_arrow_fancy_orienteering)
}

# Create different zoom levels
plot_25pct <- create_zoomed_plot(0.25)  # Middle 25%
plot_50pct <- create_zoomed_plot(0.5)   # Middle 50%
plot_75pct <- create_zoomed_plot(0.75)  # Middle 75%





grid_clipped <- st_read("grids/grid_clipped_2point5km.gpkg")
scotland <- st_read("gridExtractions/scotlandBoundary/scotland_boundary.shp")
scotland <- st_transform(scotland, crs = 27700)

###### --- grid map, 2.5km pixels --- ######
bbox <- st_bbox(scotland)
grid_clipped <- st_crop(grid_clipped, bbox)

# Find the center of the grid
grid_center <- st_centroid(st_union(grid_clipped))

# Create a buffer around the center (adjust distance as needed)
zoom_area <- st_buffer(grid_center, dist = 25000)  # 25km buffer

# Crop grid to the zoom area
grid_zoomed <- st_intersection(grid_clipped, zoom_area)

# Plot the zoomed area
ggplot() +
  geom_sf(data = grid_zoomed, linewidth = 0.1) +
  geom_sf(data = scotland, fill = NA, color = "black", linewidth = 0.3) +
  theme_bw() +
  coord_sf(expand = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.8) +
  annotation_north_arrow(location = "br", which_north = "true",
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"), width = unit(1, "cm"))


# Calculate the actual center of the grid
grid_center <- st_centroid(st_union(grid_clipped))
center_coords <- st_coordinates(grid_center)

# Define zoom area size (adjust these values as needed)
zoom_width <- 10000   # 50km wide
zoom_height <- 10000  # 50km high

# Calculate bounds around center
zoom_bbox <- c(
  xmin = center_coords[1] - zoom_width/2,
  xmax = center_coords[1] + zoom_width/2,
  ymin = center_coords[2] - zoom_height/2,
  ymax = center_coords[2] + zoom_height/2
)

zoomed_plot <- ggplot() +
  geom_sf(data = grid_clipped, linewidth = 0.1) +
  geom_sf(data = scotland, fill = NA, color = "black", linewidth = 0.3) +
  theme_bw() +
  coord_sf(
    xlim = c(zoom_bbox["xmin"], zoom_bbox["xmax"]),
    ylim = c(zoom_bbox["ymin"], zoom_bbox["ymax"]),
    expand = FALSE
  ) +
  annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.8) +
  annotation_north_arrow(location = "br", which_north = "true",
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"), width = unit(1, "cm"))

ggsave(filename = "grid_plot_center_zoom.png", 
       plot = zoomed_plot, width = 6, height = 6, 
       units = "in", dpi = 300, bg = "transparent")



### zoom to specific coordinates
bbox <- st_bbox(scotland)
grid_clipped <- st_crop(grid_clipped, bbox)

# Your coordinates in decimal degrees
manual_zoom_dd <- c(
  xmin = -4.9698,  # West
  xmax = -4.6466,  # East
  ymin = 56.8831,  # South
  ymax = 57.1641   # North
)
zoom_bbox_wgs84 <- st_bbox(manual_zoom_dd, crs = 4326)
zoom_bbox_sf <- st_as_sfc(zoom_bbox_wgs84)

# Transform to your project CRS (27700)
zoom_bbox_27700 <- st_transform(zoom_bbox_sf, crs = 27700)
zoom_bbox_coords <- st_bbox(zoom_bbox_27700)

# Use the transformed coordinates
library(ggplot2)
library(ggspatial)
library(sf)

# Only the zoomed portion of grid
zoomed_cells <- st_crop(grid_clipped, zoom_bbox_27700)
# Crop grid to the zoomed bounding box
grid_zoomed <- st_crop(grid_clipped, zoom_bbox_27700)

# Optional: crop Scotland boundary as well
scotland_zoomed <- st_crop(scotland, zoom_bbox_27700)

# Plot
grid_zoomed <- st_crop(grid_clipped, zoom_bbox_27700)
scotland_zoomed <- st_crop(scotland, zoom_bbox_27700)

zoomed_plot_27700 <- ggplot() +
  geom_sf(data = grid_zoomed, linewidth = 0.5, fill = NA, color = "black") +
  geom_sf(data = scotland_zoomed, fill = NA, color = "black", linewidth = 0.5) +
  coord_sf(
    xlim = c(zoom_bbox_coords["xmin"], zoom_bbox_coords["xmax"]),
    ylim = c(zoom_bbox_coords["ymin"], zoom_bbox_coords["ymax"]),
    expand = FALSE
  ) +
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "transparent", colour = NA),
    plot.background  = element_rect(fill = "transparent", colour = NA),
    panel.grid.major = element_blank(),  # remove major grid lines
    panel.grid.minor = element_blank()   # remove minor grid lines
  ) +
  annotation_scale(location = "bl", width_hint = 0.3, text_cex = 0.8) +
  annotation_north_arrow(location = "br", which_north = "true",
                         style = north_arrow_fancy_orienteering,
                         height = unit(1, "cm"), width = unit(1, "cm"))

zoomed_plot_27700




# Save as transparent PNG
ggsave(
  filename = "zoomed_grid_27700.png",
  plot = zoomed_plot_27700,
  width = 6,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "transparent"
)















# ======================================================================

# If you have many dataframes, use a list approach
dataframes_list <- list(survey_sf, cs_sf)  # Add more if needed
names_list <- c("survey", "cs")  # Names for source column

combined_data <- map2_dfr(dataframes_list, names_list, function(df, name) {
  st_join(df, grid_clipped, join = st_within) %>%
    st_drop_geometry() %>%
    mutate(source = name)
})

points_with_grid <- grid_clipped %>%
  left_join(combined_data, by = "grid_id")

ggplot() +
  geom_sf(data = grid_clipped, linewidth = 0.1) +
  geom_sf(data = points_with_grid, size = 2) +
  theme_minimal()

ggplot() +
  geom_sf(data = grid_clipped, linewidth = 0.1) +
  geom_sf(data = points_with_grid %>% 
            filter(!is.na(source)) %>% 
            st_as_sf(coords = c("geometry"), crs = st_crs(grid_clipped)),
          aes(color = source), 
          size = 3, 
          alpha = 0.7,
          show.legend = "point") +
  
  # # add grid ID labels
  # geom_sf_text(data = grid_clipped, 
  #              aes(label = grid_id), 
  #              size = 0.5, 
  #              color = "darkblue",
  #              fontface = "bold") +
  
  scale_color_manual(values = c("survey" = "red", "cs" = "blue")) +
  
  labs(       color = "Data Source",
       caption = paste("Total points:", nrow(combined_data), 
                       "| Survey:", sum(combined_data$source == "survey"),
                       "| Citizen Science:", sum(combined_data$source == "cs"))) +
  
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 12),
    legend.position = "bottom",
    panel.grid = element_blank()
  )

# title = "Spatial Join of Multiple Datasets with Grid",
# subtitle = "Points from two datasets mapped to grid cells"

# Count points per grid cell by source
summary_stats <- combined_data %>%
  group_by(grid_id, source) %>%
  summarise(point_count = n(), .groups = "drop") %>%
  complete(grid_id = 1:max(grid_id), source, fill = list(point_count = 0))

# simplified versions with only geometry for mapping
survey_geom <- st_as_sf(data.frame(geometry = st_geometry(survey_sf)))
cs_geom <- st_as_sf(data.frame(geometry = st_geometry(cs_sf)))

mapview(grid_clipped, color = "blue", layer.name = "Grid Lines") +
  mapview(survey_geom, col.regions = "red", layer.name = "survey") +
  mapview(cs_geom, col.regions = "green", layer.name = "citizen science")


## --- SHIFT THE GRID ALIGNMENT ---
#> handle points that fall on grid cell edges or boundaries
# visualize where points fall relative to cell centers
# Compute centroids of grid cells
grid_centroids <- st_centroid(grid_clipped)


combined_data <- map2_dfr(dataframes_list, names_list, function(df, name) {
  st_join(df, grid_clipped, join = st_within) %>%
    # st_drop_geometry() %>%
    mutate(source = name)
})

library(sf)
library(dplyr)

### MEAN OFFSET ###
# Get centroids as a data.frame with coords
grid_centroids <- st_centroid(grid_clipped) %>%
  mutate(x_c = st_coordinates(.)[,1],
         y_c = st_coordinates(.)[,2]) %>%
  st_drop_geometry() %>%
  select(grid_id, x_c, y_c)

# Attach centroid coordinates to points
points_with_centroid <- combined_data %>%
  left_join(grid_centroids, by = "grid_id")

# Get point coordinates
pt_coords <- st_coordinates(points_with_centroid)

# Compute offsets: (point - cell centroid)
offsets <- cbind(pt_coords[,1] - points_with_centroid$x_c,
                 pt_coords[,2] - points_with_centroid$y_c)

# Average offset (the shift you need)
mean_offset <- colMeans(offsets)
mean_offset
# e.g. points tend to sit 74m west and 10m south of their cell centers

grid_shifted <- st_geometry(grid_clipped) + mean_offset
grid_shifted <- st_set_geometry(grid_clipped, grid_shifted)

st_crs(grid_shifted) <- st_crs(grid_clipped)


### OPTIMAL SHIFT ###
evaluate_shift <- function(points, grid, dx, dy) {
  # Shift the grid
  shifted <- st_geometry(grid) + c(dx, dy)
  shifted <- st_set_geometry(grid, shifted)
  
  # Compute centroids of shifted grid
  centroids <- st_centroid(shifted) %>%
    mutate(x_c = st_coordinates(.)[,1],
           y_c = st_coordinates(.)[,2]) %>%
    st_drop_geometry() %>%
    select(grid_id, x_c, y_c)
  
  # Join centroids to points
  pts <- points %>%
    left_join(centroids, by = "grid_id")
  
  # Distances from points to centroids
  pt_coords <- st_coordinates(pts)
  dists <- sqrt((pt_coords[,1] - pts$x_c)^2 + (pt_coords[,2] - pts$y_c)^2)
  
  # Return mean distance
  mean(dists, na.rm = TRUE)
}

# Candidate shifts
candidates <- expand.grid(dx = seq(-500, 500, 100),
                          dy = seq(-500, 500, 100))

# Evaluate each candidate
candidates$mean_dist <- mapply(
  function(dx, dy) evaluate_shift(combined_data, grid_clipped, dx, dy),
  candidates$dx, candidates$dy
)

best <- candidates[which.min(candidates$mean_dist), ]
best

grid_best <- st_geometry(grid_clipped) + c(best$dx, best$dy)
grid_best <- st_set_geometry(grid_clipped, grid_best)

# Plot region / original grid
plot(st_geometry(grid_clipped), border = "grey70", col = NA, main = "Points with Original vs Shifted Grid")

# Add points
plot(st_geometry(combined_data), col = "red", pch = 20, cex = 0.7, add = TRUE)

# Add best-shifted grid (different color)
plot(st_geometry(grid_best), border = "blue", lwd = 1.5, col = NA, add = TRUE)

# Legend
legend("topright", legend = c("Points", "Original grid", "Shifted grid"),
       col = c("red", "grey70", "blue"), pch = c(20, NA, NA),
       lty = c(NA, 1, 1), lwd = c(NA, 1, 1), bty = "n")

grid_clipped_translated <- grid_best
st_crs(grid_clipped_translated) <- st_crs(grid_clipped)

mapview(grid_shifted, color = "purple", layer.name = "Shifted (mean) Grid Lines") +
  mapview(grid_clipped_translated, color = "pink", layer.name = "Shifted (optimal) Grid Lines") +
  mapview(survey_geom, col.regions = "red", layer.name = "survey") +
  mapview(cs_geom, col.regions = "green", layer.name = "citizen science") +
  mapview(grid_clipped, color = "blue", layer.name = "Old Grid Lines") 
  
  
  
### CREATE NEW 1KM GRID CENTERED OVER POINTS
#> generate a new 1 km grid, with the grid origin (anchor point) chosen so
#>  that points fall as close as possible to their respective cell centers.
#>  By default, the anchor is (0,0) in the CRS, change the anchor 
#>  (e.g., shift by 100 m east, 200 m north), you get a different grid, 
#>  but still 1 km cells.
#>  workflow: search over candidate anchors → build grids → evaluate how well
#>   points align with cell centers → pick the best.
#>   

# Function to make a clipped grid for a given anchor:
make_clipped_grid <- function(region, cellsize = 1000, offset = c(0,0)) {
  st_make_grid(
    region,
    cellsize = cellsize,
    offset = offset,       # where to start the grid
    what = "polygons",
    square = TRUE
  ) |>
    st_as_sf() |>
    st_intersection(region) |>
    mutate(grid_id = row_number())
}

# Function to evaluate alignment (mean point–centroid distance)
evaluate_grid <- function(points, bb, dx, dy, cellsize = 1000) {
  # point coordinates
  pt_coords <- st_coordinates(points)
  
  # shifted grid origin
  x0 <- bb["xmin"] + dx
  y0 <- bb["ymin"] + dy
  
  # compute the centroid of the cell each point falls in
  cell_x <- floor((pt_coords[,1] - x0) / cellsize) * cellsize + x0 + cellsize/2
  cell_y <- floor((pt_coords[,2] - y0) / cellsize) * cellsize + y0 + cellsize/2
  
  # distances to centroids
  dists <- sqrt((pt_coords[,1] - cell_x)^2 + (pt_coords[,2] - cell_y)^2)
  
  mean(dists, na.rm = TRUE)
}

bb <- st_bbox(scotland)

candidates <- expand.grid(dx = seq(-500, 500, 100),
                          dy = seq(-500, 500, 100))

candidates$mean_dist <- mapply(
  function(dx, dy) evaluate_grid(combined_data, bb, dx, dy, cellsize = 1000),
  candidates$dx, candidates$dy
)


# Find the best offset
best <- candidates[which.min(candidates$mean_dist), ]
best

# Build final optimized grid
best_grid <- make_clipped_grid(
  scotland, 1000,
  offset = c(bb["xmin"] + best$dx, bb["ymin"] + best$dy)
)

plot(st_geometry(region_shape), border = "black", col = NA,
     main = "Optimized 1km Grid with Points")
plot(st_geometry(best_grid), border = "blue")
plot(st_geometry(points), col = "red", pch = 20, cex = 0.7, add = TRUE)

st_write(best_grid, "grids/best_grid.shp")

mapview(best_grid, color = "purple", layer.name = "New Grid Lines") +
  mapview(survey_geom, col.regions = "red", layer.name = "survey") +
  mapview(cs_geom, col.regions = "green", layer.name = "citizen science") +
  mapview(grid_clipped, color = "blue", layer.name = "Old Grid Lines") 



