#####################################################
# Prepare cost surfaces -----------------------------------------
#####################################################

# Libraries
library(sf)
library(raster)
library(mapview)
library(fasterize)
library(doParallel)
doParallel::registerDoParallel(parallel::detectCores() - 1)


# Load spatial vectors
load("data/spatialVectors.rda")
land_ca <- raster::raster("data/landcover_ca_30m.tif")



#####################################################################
# Define cost parameters
#####################################################################

# number of ARUs to be deployed for each hexagon
params <- list(
  nb_arus = 5.0,
  # Truck cost
  truck_buffer = 1000, # meters
  truck_cost_per_day = 600.0,
  truck_n_crews	= 2.0,
  truck_arus_per_crew_per_day	= 5.0,
  # ATV cost
  atv_buffer = 1000,
  atv_cost_per_day = 1200.0,
  atv_n_crews	= 2.0,
  atv_arus_per_crew_per_day	= 3.0,
  # Helicopter
  helicopter_cost_per_hour = 1250.0,
  helicopter_max_km_from_base	= 150.0,
  helicopter_base_setup_cost_per_km =	9.0,
  helicopter_l_per_hour = 160.0,
  helicopter_crew_size = 4.0,
  helicopter_aru_per_person_per_day	= 5.0,
  helicopter_relocation_speed = 180.0,
  helicopter_airport_cost_per_l	= 1.3,
  helicopter_base_cost_per_l	= 5.0,
  helicopter_2nd_base_cost_per_l = 10.0,
  helicopter_hours_flying_within_sa_per_day	= 5.0
)


#####################################################################
# Roads (truck)
#####################################################################

# Buffer around all roads
road_buff <- sf::st_buffer(roads, dist = params$truck_buffer)

# Standardize all polygons type to polygon (fasterize accepts only ONE type of geometry)
road_buff <- sf::st_cast(road_buff, 'POLYGON')

# rasterize road buffer
road_buff_r <- fasterize::fasterize(road_buff, land_ca)

# Add a fixed cost depending on parameters for the buffer
truckFixedCost <- params$truck_cost_per_day * params$nb_arus / params$truck_arus_per_crew_per_day
road_cost <- raster::reclassify(road_buff_r, cbind(1, truckFixedCost))

# clean memory
rm(road_buff, road_buff_r, roads)



#####################################################################
# ATV (trails)
#####################################################################

# Buffer around all ATV trails
atv_buff <- sf::st_buffer(trails, dist = params$atv_buffer)

# Standardize all polygons type to polygon (fasterize accepts only ONE type of geometry)
atv_buff <- sf::st_cast(atv_buff, 'POLYGON')

# rasterize road buffer
atv_buff_r <- fasterize::fasterize(atv_buff, land_ca)

# Add a fixed cost for the buffer
atvFixedCost <- params$atv_cost_per_day * params$nb_arus / params$atv_arus_per_crew_per_day
atv_cost <- raster::reclassify(atv_buff_r, cbind(1, atvFixedCost))

# clean memory
rm(atv_buff, atv_buff_r, trails)



#####################################################################
# Helicopter access
#####################################################################

# Get centroid of hexagons
# Helicoter cost is calculated once for the whole hexagon based on its centroid (not pixel based)
hexa_hel <- hexa
hexa_hel$geometry <- hexa %>%
                      sf::st_centroid() %>%
                      sf::st_geometry() %>%
                      sf::st_transform(sf::st_crs(aeroports))

# Train rail is used as a source of fuel (no need to double travel and basecamp, so they are considered as aerports)
# Get train rail in which is inside buffer of X km around aeroport (X is max distance helicopter can fly)
aeroports_bf <- sf::st_buffer(aeroports, dist = params$helicopter_max_km_from_base * 1000)

# select all hexagons in which cetroid is inside this buffer
train_closeAeroport <- sf::st_intersection(train, aeroports_bf)

# Calculate matrix of distance between centroid and aeroports
mat_dist <- sf::st_distance(x = hexa_hel, y = aeroports)
hexa_hel$aeroportDist <- apply(mat_dist, 1, min)

# Same as above but for the train rail
mat_dist <- sf::st_distance(x = hexa_hel, y = train_closeAeroport)
hexa_hel$trainDist <- apply(mat_dist, 1, min)

# minimum distance between aeroports and train
hexa_hel$aeroportDist <- pmin(hexa_hel$aeroportDist, hexa_hel$trainDist)
hexa_hel <- hexa_hel[, !(names(hexa_hel) %in% 'trainDist')]


# Add cost depending on distance (meters)
heliCost <- function(distance, pars)
{
  # Transform distance to km
  distance <- distance/1000

  # Cost per litre of fuel (depending on distance)
  Cl <- ifelse(distance < pars$helicopter_max_km_from_base,
                pars$helicopter_airport_cost_per_l,
                ifelse(distance < (2 * pars$helicopter_max_km_from_base),
                  pars$helicopter_base_cost_per_l,
                  pars$helicopter_2nd_base_cost_per_l))

  # Cost per hour of helicopter (rent + fuel)/h
  Ch <- pars$helicopter_l_per_hour * Cl + pars$helicopter_cost_per_hour

  # Cost per Km of helicopter (rent + fuel)
  Cd <- Ch/pars$helicopter_relocation_speed

  # Hours per ARUs deployed (depending on h/day of helicopter, number of crew, and AURs per crew)
  H_arus <- pars$helicopter_hours_flying_within_sa_per_day/(pars$helicopter_crew_size * pars$helicopter_aru_per_person_per_day)

  # If basecamp needed, calculate:
  # - distance from (A)report to (B)asecamp
  # - distance from (B)asecamp to (H)exagon
  # - Cb - cost of flying from aerport to basecamp
  if(distance >= pars$helicopter_max_km_from_base)
  {
    distanceBH <- distance - pars$helicopter_max_km_from_base
    distanceAB <- distance - distanceBH
    
    # Cost of flying from Aeroport to basecamp
    Cb <- (distanceAB * pars$helicopter_base_setup_cost_per_km) + 
          (2 * distanceAB/pars$helicopter_relocation_speed * pars$helicopter_l_per_hour * Cl)
  }else{
    Cb <- 0
  }

  # Cost within hexagon (number of ARUs * Hours per ARU * cost per hour)
  Csa <- pars$nb_arus * H_arus * Ch

  # Cost of flying to the hexagon (either from Aerport OR basecamp)
  if(distance >= pars$helicopter_max_km_from_base)
  {
    Cf <- 2 * distanceBH * Cd
  }else{
    Cf <- 2 * distance * Cd
  }

  # final cost to get to a hexagon
  Cost_heli <- Csa + Cf + Cb
}

hexa_hel$heliCost <- sapply(hexa_hel$aeroportDist, heliCost, params)

# clean memory
rm(mat_dist)



#####################################################################
# Stack the two layers cost (roads + trails)
#####################################################################

# stack all three layers
cost_stack <- raster::stack(road_cost, atv_cost)

# get minimum cost for each cell
cost_min <- raster::calc(cost_stack, min)

# save
raster::writeRaster(cost_min, filename = 'data/min_cost.tif', format = 'GTiff')
saveRDS(hexa_hel, 'data/hexa_hel_cost.RDS')





#####################################################################
# Extract mean cost for each hexagon
#####################################################################

# load min cost raster calculated above
cost_min <- raster::raster('data/min_cost.tif')
hexa_hel <- readRDS('data/hexa_hel_cost.RDS')

for (id in unique(districts$ECOREGION))
{
    cat('Ecoregion', id, '-- Process', which(unique(districts$ECOREGION) == id), 'on', length(unique(districts$ECOREGION)), '\n')
    
    # Select ecoregion
    ecoregion <- districts[districts$ECOREGION == id, ]

    # Load hexagons of ecoregion (generated from `script 02_probHab.R`)
    hexa_ecoregion <- sf::st_read(paste0('output/', tolower(gsub(' ', '_', unique(ecoregion$REGION_NAM))), '_', id, '/hexa_phab.shp'), quiet = TRUE)

    # Crop cost raster to ecoregion
    print('Start to calculate the cost-based inclusion probabilities')
    costMin_ecoregion <- mask(crop(cost_min, as(hexa_ecoregion, "Spatial")), as(hexa_ecoregion, "Spatial"))
    
    # Extract habitat values within each hexagons -- FOR CA
    costSum <- foreach(i = seq_len(nrow(hexa_ecoregion)), .packages = "raster") %dopar% {    
      hexValues <- raster::extract(costMin_ecoregion, hexa_ecoregion[i, ])
      costTab <- table(hexValues[[1]])
      propTransport <- setNames(rep(0, 2), c(as.character(params$truck_cost_per_day),
                                           as.character(params$atv_cost_per_day)))
      if(nrow(costTab) > 0)
        propTransport[names(costTab)] <- costTab/length(hexValues[[1]])
      
      # add helicopter (vector name is the helicopter cost for the hexagon to facilitate calculation)
      heli_cost_hex_i <- as.character(subset(hexa_hel, ET_ID == hexa_ecoregion[i, ]$ET_ID)$heliCost)
      propTransport[heli_cost_hex_i] <- 1 - sum(propTransport)

      # calculate weighted cost
      sum(as.numeric(names(propTransport)) * propTransport)
    }

    # Assigning new columns with hab prob
    hexa_ecoregion$costSum <- unlist(costSum)

    # Calculate the root squared of inversed cost to get a probabilistic value
    inv_sqrt <- 1 / (sqrt(unlist(costSum)))
    hexa_ecoregion$costProb <- inv_sqrt / sum(inv_sqrt)
    
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
