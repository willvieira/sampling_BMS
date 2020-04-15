#####################################################
# Prepare cost surfaces -----------------------------------------
#####################################################

# Libraries
library(sf)
library(mapview)
library(raster)
library(doParallel)
library(tictoc)

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
registerDoParallel(3)
# Compute prevalence for quebec landcover
tic("roads_area", log = TRUE)
val <- foreach(i = seq_len(length(ref_grid)), .combine=c) %dopar% { min(st_distance(ref_grid[i], roads_area)) }
toc()
df_roads <- data.frame(st_coordinates(ref_grid), val=val)
rs_cost_roads <- rasterFromXYZ(df_roads, crs = st_crs(area)$proj4string)

# Compute euclidean distance from aeroports
tic("aeroports_area", log = TRUE)
val <- foreach(i = seq_len(length(ref_grid)), .combine=c) %dopar% { min(st_distance(ref_grid[i], aeroports_area))}
toc()
df_aeroports <- data.frame(st_coordinates(ref_grid), val=val)
rs_aeroports <- rasterFromXYZ(df_aeroports, crs = st_crs(area)$proj4string)

# Compute euclidean distance from trails
tic("trails_area", log = TRUE)
val <- foreach(i = seq_len(length(ref_grid)), .combine=c) %dopar% { min(st_distance(ref_grid[i], trails_area))}
toc()
df_trails <- data.frame(st_coordinates(ref_grid), val=val)
rs_trails <- rasterFromXYZ(df_trails, crs = st_crs(area)$proj4string)

cost_stack <- stack(rs_roads, rs_aeroports, rs_trails)

# Compute cost raster using density function
# Convert to a class exploitable by maptools and spatstat
# psp_roads <-  as.psp(as(roads_area, "Spatial"))
# psp_trails <-  as.psp(as(trails_area, "Spatial"))
# ppp_aeroport <-  as.ppp(as(aeroports_area, "Spatial"))

# cost_roads <- density(psp_roads, sigma = 1e3)
# cost_trails <- density(psp_trails, sigma = 40e3)
# cost_aeroports <- density(ppp_aeroport, sigma = 100e3)

# Coerce density back to raster
rs_cost_roads <- rasterFromXYZ(as.data.frame(cost_roads), crs = st_crs(roads_area)$proj4string)
rs_cost_trails <- rasterFromXYZ(as.data.frame(cost_trails), crs = st_crs(trails_area)$proj4string)
rs_cost_aeroports <- rasterFromXYZ(as.data.frame(cost_aeroports), crs = st_crs(aeroports_area)$proj4string)

#####################################################
# Prepare cost surfaces -----------------------------------------
#####################################################
