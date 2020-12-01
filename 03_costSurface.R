#####################################################
# Prepare cost surfaces -----------------------------------------
#####################################################

# Libraries
library(sf)
library(mapview)
library(raster)
library(fasterize)
library(doParallel)
doParallel::registerDoParallel(parallel::detectCores() - 2)


# Load spatial vectors
load("data/spatialVectors.rda")
land_ca <- raster::raster("data/landcover_ca_30m.tif")


# crop lines (roads, trails, aeroports) for the study area
roads_area <- sf::st_crop(roads, area)
trails_area <- sf::st_crop(trails, area)
aeroports_area <- sf::st_crop(aeroports, area)

# View spatial cost objects
mapview(roads_area) + mapview(trails_area) + mapview(aeroports_area)



#####################################################################
# Roads
#####################################################################

# 2km buffer around all roads
road_buff <- sf::st_buffer(roads_area, dist = 2000)

# rasterize road buffer
road_buff_r <- fasterize::fasterize(road_buff, land_ca)

# Add a fixed cost of 100 $ for the 2km buffer
road_cost <- raster::reclassify(road_buff_r, cbind(1, 100))

# clean memory
rm(road_buff, road_buff_r, roads, roads_area)



#####################################################################
# Trails
#####################################################################

# 300 m buffer around all roads
trail_buff <- sf::st_buffer(trails_area, dist = 300)

# rasterize road buffer
trail_buff_r <- fasterize::fasterize(trail_buff, land_ca)

# Add a fixed cost of 200 $ for the 300 m buffer
trail_cost <- raster::reclassify(trail_buff_r, cbind(1, 200))

# clean memory
rm(trail_buff, trail_buff_r, trails, trails_area)



#####################################################################
# Helicopter access
#####################################################################

heli_buff <- trail_cost
values(heli_buff) <- NA

# Two classes of distance (45, 150)
# Everything farther to 250 km will have probability of zero (not yet)
heli_dist <- raster::distanceFromPoints(heli_buff, as(aeroports_area[['Shape']], 'Spatial'))

# Add fixed cost depending on distance classes
heli_cost <- raster::calc(heli_dist, function(x) {
  ifelse(x <= 45000, 1000, ifelse(x < 150000, ((x / 1000) * 2 * 16.23) - 460.70, ((x / 1000) * 2 * 16.23) + 239.04))
})

# clean memory
rm(heli_dist, heli_buff)



#####################################################################
# Stack the three layers cost (roads, trails, and helicopters)
#####################################################################

# stack all three layers
cost_stack <- raster::stack(road_cost, trail_cost, heli_cost)

# get minimum cost for each cell
cost_min <- raster::calc(cost_stack, min)

# save
raster::writeRaster(cost_min, filename = 'data/min_cost.tif', format = 'GTiff')
rm(cost_stack)





#####################################################################
# Extract mean cost for each exagon
#####################################################################

# load min cost raster calculated above
cost_min <- raster::raster('data/min_cost.tif')


for (id in unique(districts$ECOREGION)) {

    cat('Ecoregion', id, '-- Process', which(unique(districts$ECOREGION) == id), 'on', length(unique(districts$ECOREGION)), '\n')
    
    # Select ecoregion
    ecoregion <- districts[districts$ECOREGION == id, ]

    # Load hexagons of ecoregion (generated from `script 02_probHab.R`)
    hexa_ecoregion <- sf::st_read(paste0('output/', tolower(gsub(' ', '_', unique(ecoregion$REGION_NAM))), '_', id, '/hexa_phab.shp'))

    # Extract habitat values within each hexagons -- FOR CA
    print('Start to extract habitat values for CA')
    costMin <- foreach(i = seq_len(nrow(hexa_ecoregion)), .packages = "raster", .combine = c) %dopar% {    
        raster::extract(cost_min, hexa_ecoregion[i, ], fun = mean)
    }

    # Assigning new columns with hab prob
    hexa_ecoregion$cost_min <- costMin
    #hexa_ecoregion$hab_qc <- hab_qc

    # Add new phab columns in attributes table
    print('Saving')

    # Writing results in corresponding folder
    write_sf(hexa_ecoregion, paste0('output/', tolower(gsub(' ', '_', unique(ecoregion$REGION_NAM))), '_', id, '/hexa_cost.shp'))
}

stopImplicitCluster()




#####################################################################
# Inclusion probability based on distance only (Not used)
#####################################################################

# # Make reference grid (resolution: 1K)
# # Protocol suggest to stick to lcc/modis resolution but require lot of RAM
# ref_grid <- st_make_grid(area, 10000, what="centers")

# # Compute euclidean distance from roads
# # Open Cluster
# registerDoParallel(detectCores()-2)

# # Compute euclidean distance from roads
# val <- foreach(i = seq_len(length(ref_grid)), .combine=c) %dopar% { min(st_distance(ref_grid[i], roads_area)) }
# df_roads <- data.frame(st_coordinates(ref_grid), val=val)
# rs_roads <- rasterFromXYZ(df_roads, crs = st_crs(area)$proj4string)

# # Compute euclidean distance from aeroports
# val <- foreach(i = seq_len(length(ref_grid)), .combine=c) %dopar% { min(st_distance(ref_grid[i], aeroports_area))}
# df_aeroports <- data.frame(st_coordinates(ref_grid), val=val)
# rs_aeroports <- rasterFromXYZ(df_aeroports, crs = st_crs(area)$proj4string)

# # Compute euclidean distance from trails
# val <- foreach(i = seq_len(length(ref_grid)), .combine=c) %dopar% { min(st_distance(ref_grid[i], trails_area))}
# df_trails <- data.frame(st_coordinates(ref_grid), val=val)
# rs_trails <- rasterFromXYZ(df_trails, crs = st_crs(area)$proj4string)


# # Get minimum cost for each point
# cost_stack <- raster::stack(rs_roads, rs_aeroports, rs_trails)
# min_cost <- raster::calc(cost_stack, min)

# # obtain inverse square root of cost
# invsqr_cost <- calc(min_cost, function(x) {
#   1 / (sqrt(x))
# }) # invsqrcost
