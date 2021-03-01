### Title: Preparing spatial layers for analysis
### Author: Steve Vissault
### Last edited by Will Vieira on feb 21, 2020


# Import Packages ---------------------------------------------------------
library(raster)
library(sf)
library(mapview)

# Parallel computing 
raster::beginCluster()


# Set raster tmp folder
# rasterOptions(tmpdir = "/data/sviss/tmp")

###############################################
# Prepare spatial layers ---------------------------------------------------------
###############################################

######## IMPORT: study area

# Read study area
area <- read_sf("rawLayers/SOBQ_ATLAS_PERI_Micro_W.shp")

# Read ecoregions
districts <- read_sf("rawLayers/Ecoregions_Canada_20201209_W.shp")
#mapview::mapview(districts, zcol = 'ECOREGION')

# Read Hexagons
# 414163 unique ET_ID
hexa <- read_sf("rawLayers/HexagonesCND_SOBQ_W2.shp")

######## IMPORT: landcover Canada and Canada adapted from land Quebec
# Read LCC2015: Habitat proportion
land_ca <- raster("rawLayers/CAN_LC_2015_CAL.tif")

# Read LCC2015 adapted from land_Qc (combined info, needs reclassification using `rawLayers/ClassesVegSOBQ_combineV1Fev2021.csv`)
land_caqc <- raster("rawLayers/Landcover_Combine.tif")

# Layer defining the limits of land_caqc within Quebec
limits_caqc <- raster('rawLayers/UT2017_Extent.tif')


######## IMPORT: layers for cost analysis

# Read roads
roads_qc <- st_read("rawLayers/Reseau_routier_SOBQ_ATLAS_W.shp")
roads_lb <- st_read("rawLayers/Reseau_routier_Labrador_50kW.shp")

# Read trails 
trails <- st_read("rawLayers/Motoneiges_Qc_W.shp")

# Read Aeroports
aeroports_qc <- st_read("rawLayers/Aeroports_Qc_20KW.shp")
aeroports_lb <- st_read("rawLayers/Aeroport_Labrador_50KW.shp")

# Read train
train <- st_read("rawLayers/Trains_TshiuetinQCLAB_W.shp")


####### TRANSFORM: reproject on land_ca

# Reproject all spatial layers under land_ca cover projection
area <- st_transform(area, st_crs(land_ca))
districts <- st_transform(districts, st_crs(land_ca))
hexa <- st_transform(hexa, st_crs(land_ca))
roads_qc <- st_transform(roads_qc, st_crs(land_ca))
roads_lb <- st_transform(roads_lb, st_crs(land_ca))
trails <- st_transform(trails, st_crs(land_ca))
aeroports_qc <- st_transform(aeroports_qc, st_crs(land_ca))
aeroports_lb <- st_transform(aeroports_lb, st_crs(land_ca))
train <- st_transform(train, st_crs(land_ca))


####### TRANSFORM: Crop and group on study area

# Filter hexagons for study area (rule: hexagons must have centroid inside study area)
hexa_cent <- hexa
hexa_cent$geometry <- hexa %>% sf::st_centroid() %>% sf::st_geometry()
hexa_cent <- hexa_cent[which(st_intersects(hexa_cent, area, sparse = FALSE)), ]
hexa <- hexa[which(hexa$ET_Index %in% hexa_cent$ET_Index), ]

districts <- st_intersection(districts, area)
# Remove specific ecoregions
ecoregionsToRm <- c(97, 118, 82, 104)
districts <- districts[which(!districts$ECOREGION %in% ecoregionsToRm), ]



# We are not croping aeroports and train as some elements of these layers that
# are outside study area, but close, may also be used to access sites
# For roads and trails, we will expand (buffer) study area by 10 km to get roads and trails outside but close to the study area
area_expanded <- st_buffer(area, dist = 10000)
roads_qc <- st_crop(roads_qc, area_expanded)
roads_lb <- st_crop(roads_lb, area_expanded)
trails <- st_crop(trails, area_expanded)


####### TRANSFORM: filter aeroports
## Infrastructure: either 'Aéroport', 'Héliport' or 'Aérodrome'
## Carburant == 'OUI'
aeroports_qc <- aeroports_qc[which(aeroports_qc$typeinfras %in% c('Aéroport', 'Héliport', 'Aérodrome')), ]
# aeroports to keep whatever if carburant is Y or N
aeroToKeep <- aeroports_qc[which(aeroports_qc$nomexploit %in% c('Hydro-Québec', 'Administration Régionale Kativik')), ]
# Get all aeroports with carburant
aeroports_qc <- aeroports_qc[which(aeroports_qc$carburant == 'OUI'), ]
# add the remaining aerports to keep
aeroports_qc <- rbind(aeroports_qc, aeroToKeep)


####### TRANSFORM: Group Quebec and Labrador layers
# Aeroports (keep only ID in attribute table, and add qc OR lb to id of each aeroport)
aeroports_qc <- setNames(aeroports_qc['idobj'], c('OBJECTID', 'geometry')) # select only 'idobj' column and change its name to OBJECTID
aeroports_qc$OBJECTID <- paste0(aeroports_qc$OBJECTID, '_qc')
aeroports_lb <- aeroports_lb['OBJECTID']
aeroports_lb$OBJECTID <- paste0(aeroports_lb$OBJECTID, '_lb')
aeroports <- rbind(aeroports_qc, aeroports_lb)

# roads
roads_qc <- roads_qc['OBJECTID']
roads_qc$OBJECTID <- paste0(roads_qc$OBJECTID, '_qc')
roads_lb <- roads_lb['OBJECTID']
roads_lb$OBJECTID <- paste0(roads_lb$OBJECTID, '_lb')
roads <- rbind(roads_qc, roads_lb)


###############################################
# Crop and reclassify landcover (CA & CAQC) ---------------------------------------------------------
###############################################

##### LCC CANADA #####

# Crop land_ca with area (reduce memory alloc)
area_lcc_ca <- st_transform(area, projection(land_ca))
land_ca <- mask(crop(land_ca, area_lcc_ca), area_lcc_ca)


# From Ontario script (see ON/03-PrepareHexagons_LCC2015_Ontario.R)
# Removing Snow and ice, water, Urban, and cropland to NA
# Reclassify
from_to <- matrix(c(0:19, c(NA, 1:14, NA, 16, rep(NA, 3))), ncol = 2)
land_ca <- reclassify(land_ca, from_to)



##### LCC CANADA adapted from QUBEC #####

# Mask land_caqc using limits_caqc (everything outside limits_caqc will be NA)
land_caqc <- overlay(land_caqc, limits_caqc, fun = function(x, y) ifelse(y < 1, NA, x))

# Load classes from combine file which combines canada and quebec information
classVegSOBQ <- read.csv("rawLayers/ClassesVegSOBQ_combineV1Fev2021.csv")
# Removing Snow and ice, water, Urban, and cropland to NA
classVegSOBQ$ClassesSOBQVEG_V1_Fev2021[which(classVegSOBQ$ClassesSOBQVEG_V1_Fev2021 %in% c(0, 15, 17, 18, 19))] <- NA
# Transform into a matrix
from_to <- as.matrix(classVegSOBQ[, 2:3])
# Reclassify
land_caqc <- reclassify(land_caqc, from_to)


###############################################
# Save clean spatial layers in Data -----`----------------------------------------------------
###############################################

# Save rasters/Landcovers
writeRaster(land_ca, filename = "data/landcover_ca_30m", format = "GTiff")
writeRaster(land_caqc, filename = "data/landcover_caqc_30m", format = "GTiff")

# Save vector data
save(area, districts, hexa, roads, trails, aeroports, train, file = "data/spatialVectors.rda")
