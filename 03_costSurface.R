#####################################################
# Prepare cost surfaces -----------------------------------------
#####################################################

# Libraries
library(sf)
library(mapview)
library(raster)
library(doParallel)

# Load spatial vectors
load("data/spatialVectors.rda")

# crop lines (roads, trails, aeroports) for the study area
roads_area <- st_crop(roads, area)
aeroports_area <- st_crop(aeroports, area)
trails_area <- st_crop(trails, area)

# View spatial cost objects
mapview(roads_area) + mapview(aeroports_area) + mapview(trails_area)

# Make reference grid (resolution: 1K)
# Protocol suggest to stick to lcc/modis resolution but require lot of RAM
ref_grid <- st_make_grid(area, 10000, what="centers")

# Compute euclidean distance from roads
# Open Cluster
registerDoParallel(detectCores()-2)

# Compute euclidean distance from roads
val <- foreach(i = seq_len(length(ref_grid)), .combine=c) %dopar% { min(st_distance(ref_grid[i], roads_area)) }
df_roads <- data.frame(st_coordinates(ref_grid), val=val)
rs_cost_roads <- rasterFromXYZ(df_roads, crs = st_crs(area)$proj4string)

# Compute euclidean distance from aeroports
val <- foreach(i = seq_len(length(ref_grid)), .combine=c) %dopar% { min(st_distance(ref_grid[i], aeroports_area))}
df_aeroports <- data.frame(st_coordinates(ref_grid), val=val)
rs_aeroports <- rasterFromXYZ(df_aeroports, crs = st_crs(area)$proj4string)

# Compute euclidean distance from trails
val <- foreach(i = seq_len(length(ref_grid)), .combine=c) %dopar% { min(st_distance(ref_grid[i], trails_area))}
df_trails <- data.frame(st_coordinates(ref_grid), val=val)
rs_trails <- rasterFromXYZ(df_trails, crs = st_crs(area)$proj4string)

cost_stack <- stack(rs_roads, rs_aeroports, rs_trails)

