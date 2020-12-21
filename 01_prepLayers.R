### Title: Preparing spatial layers for analysis
### Author: Steve Vissault

# Import Packages ---------------------------------------------------------
library(spsurvey) # GRTS
library(raster)
library(sf)
library(dplyr)
library(mapview)
library(spatstat)
library(maptools)
library(ggspatial)

# Set raster tmp folder
# rasterOptions(tmpdir = "/data/sviss/tmp")
    
###############################################
# Prepare spatial layers ---------------------------------------------------------
###############################################

######## IMPORT: study area

# Read study area
area <- read_sf("rawLayers/SOBQ_Nord.shp")

# Read eco districts
# 17 unique ECOREGIONS
# 53 unique ECODISTRICS
districts <- read_sf("rawLayers/Ecodistricts_SOBQ_Nord.shp")
#mapview::mapview(districts, zcol = 'ECODISTRIC')

# Read Hexagons
# 414163 unique ET_ID
hexa <- read_sf("rawLayers/Hexagons_w_Attributes.shp")

######## IMPORT: landcover qc & ca

# Read LCC2015: Habitat proportion
land_ca <- raster("rawLayers/CAN_LC_2015_CAL.tif")

# Read Lancover Québec
# See https://www.donneesquebec.ca/recherche/fr/dataset/vegetation-du-nord-quebecois
# Metadata: https://TinyURL.com/uyvswpb
# land_qc <- st_read("rawLayers/veg_nord_53.gdb", layer = "veg_nord")

######## IMPORT: layers for cost analysis

# Read roads
roads <- st_read("rawLayers/Couches_SamplingR.gdb", layer = "road_segment_1")

# Read trails
trails <- st_read("rawLayers/Couches_SamplingR.gdb", layer = "Sentiers")

# Read Aeroportsq
aeroports <- st_read("rawLayers/Couches_SamplingR.gdb", layer = "Aeroports_SOBQ")

####### TRANSFORM: reproject on study area

# Reproject all spatial layers under land_qc cover projection
area <- st_transform(area, st_crs(land_ca))
districts <- st_transform(districts, st_crs(land_ca))
hexa <- st_transform(hexa, st_crs(land_ca))
roads <- st_transform(roads, st_crs(land_ca))
trails <- st_transform(trails, st_crs(land_ca))
aeroports <- st_transform(aeroports, st_crs(land_ca))

###############################################
# Crop and reclassify landcover (CA & QC) ---------------------------------------------------------
###############################################

##### LCC CANADA #####

# Crop land_ca with area (reduce memory alloc)
area_lcc_ca <- st_transform(area, projection(land_ca))
land_ca <- mask(crop(land_ca, area_lcc_ca), area_lcc_ca)

# Reproject
raster::beginCluster(n = 12) # Parallel if needed
land_ca <- projectRaster(land_ca, crs=st_crs(land_qc)$proj4string, method = "ngb")
raster::endCluster()

# From Ontario script (see ON/03-PrepareHexagons_LCC2015_Ontario.R)
# Removing Snow and ice, water, Urban, and cropland to NA
# Reclassify
from_to <- matrix(c(0:19, c(NA, 1:14, NA, 16, rep(NA, 3))), ncol = 2)
land_ca <- reclassify(land_ca, from_to)

##### LCC QUEBEC #####

# # Transform land_qc to raster
# # Select habitat type column
# land_qc <- dplyr::select(land_qc, CL_CARTO)
# # Store metadata
# md_co_ter <- data.frame(gov_code = levels(land_qc$CL_CARTO), code = 1:nlevels(land_qc$CL_CARTO))

# # Coerce character to integer
# land_qc$CL_CARTO <- as.numeric(land_qc$CL_CARTO)

# # Rasterize landcover polygons QC with same resolution than land_ca
# library(fasterize)
# land_qc <- fasterize(land_qc, raster(land_qc, res = 30), field = "CL_CARTO", fun = "first")

# # Crop land_ca with study area (reduce memory alloc)
# land_qc <- crop(land_qc, area)

# # Removing Snow and ice, water, Urban, and cropland to NA
# # Medatadata classification: see object md_co_ter
# # Prepare new surface code column
# md_co_ter$new_code <- md_co_ter$code

# # Eau (code: EAU)
# md_co_ter$new_code[which(md_co_ter$gov_code == "EAU")] <- NA
# # Suface neige (Code: NE)
# md_co_ter$new_code[which(md_co_ter$gov_code == "NE")] <- NA
# # Infra humaine (Code: IH)
# md_co_ter$new_code[which(md_co_ter$gov_code == "IH")]  <- NA
# # Ligne de transport �l�c. (Code: LTE)
# md_co_ter$new_code[which(md_co_ter$gov_code == "LTE")]  <- NA

# # Run reclassification on landcover QC
# land_qc <- reclassify(land_qc, matrix(c(md_co_ter$code, md_co_ter$new_code), ncol = 2))


###############################################
# Save clean spatial layers in Data ---------------------------------------------------------
###############################################

# Save rasters/Landcovers
writeRaster(land_ca, filename = "data/landcover_ca_30m", format = "GTiff")
writeRaster(land_qc, filename = "data/landcover_qc_30m", format = "GTiff")

# Save vector data
save(area, districts, hexa, roads, trails, aeroports, file = "data/spatialVectors.rda")
