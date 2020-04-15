### Title: Compute inclusion probability based on habitat prevalence
### Author: Steve Vissault

library(raster)
library(sf)

load("data/spatialVectors.rda")
land_ca <- raster("data/landcover_ca_30m.tif")
land_qc <- raster("data/landcover_qc_30m.tif")

###############################################
# Habitat prevalence within each ecoregions ---------------------------------------------------------
###############################################
# See page 32: Protocols/BMS_Boreal Monit_LB_DRAFT_May2_2018.docx

# Function to compute the prevalence habitat (area by habitat / area district) for each ecodistricts
prevalence_hab <- function(district, landcover){

    # Subselect area
    landcover <- mask(crop(landcover, district), district)

    # Inclusion prob
    # Count pixels by habitat type
    freq_hab <- as.data.frame(table(landcover[]))
    names(freq_hab) <- c("code", "freq")

    # Compute the area based on raster res
    freq_hab$area_hab <- freq_hab$freq * (res(landcover)[1] * res(landcover)[2])
    # Compute inclusion proba
    freq_hab$incl_prob <- 1 / nrow(freq_hab) / freq_hab$area_hab

    # Add poly ID
    freq_hab$ID_poly <- district$ECODISTRIC

    # Return the dataframe
    return(freq_hab)
}


# Parallelize function with doParallel
library(doParallel)

# Open Cluster
registerDoParallel(cores = detectCores()-1)
# Compute prevalence for quebec landcover
prev_qc_by_district <- foreach(i = seq_len(nrow(districts)), .packages = "raster") %dopar% {
    try(prevalence_hab(district = districts[i,], landcover = land_qc))
}

# Compute prevalence for canada landcover
prev_ca_by_district <- foreach(i = seq_len(nrow(districts)), .packages = "raster") %dopar% {
    try(prevalence_hab(district = districts[i,], landcover = land_ca))
}

####### Save prevalence by district
prev_all_ca <- do.call(rbind, prev_ca_by_district[which(sapply(prev_ca_by_district, class) == "data.frame")])
prev_all_qc <- do.call(rbind, prev_qc_by_district[which(sapply(prev_qc_by_district, class) == "data.frame")])
# save(prev_all_qc, prev_all_ca, file="outputs/prevalence_by_district.rda")


###############################################
# Select ecoregion ---------------------------------------------------------
###############################################

load("outputs/prevalence_by_district.rda")
# future iterator
sel_district = 300

# Schefferville district
district_scheffer <- districts[districts$ECODISTRIC == sel_district,]

# Select hexagons within the area (st_insersects preserve topology)
hexa_area <- hexa[st_intersects(hexa, district_scheffer, sparse = FALSE),]

# Select prevalence for that district
prev_district_ca <- subset(prev_all_ca, ID_poly == sel_district)
prev_district_qc <- subset(prev_all_qc, ID_poly == sel_district)

# Extract habitat values within each hexagons -- FOR CA
pb_hab_hexas_qc <- foreach(i = seq_len(nrow(hexa_area)), .packages = "raster", .combine = c) %dopar% {
    count_hab <- as.data.frame(table(extract(land_qc, hexa_area[i,])))
    
    if(nrow(count_hab) > 0){
        prev_count_hab <- merge(count_hab, prev_district_qc, by.x="Var1", by.y="code" ,all.x=TRUE)
        sum(prev_count_hab$Freq * prev_count_hab$incl_prob)
    } else {
        NA
    }
}

# Extract habitat values within each hexagons -- FOR QC
pb_hab_hexas_ca <- foreach(i = seq_len(nrow(hexa_area)), .packages = "raster", .combine = c) %dopar% {
    count_hab <- as.data.frame(table(extract(land_ca, hexa_area[i,])))
    
    if(nrow(count_hab) > 0){
        prev_count_hab <- merge(count_hab, prev_district_ca, by.x="Var1", by.y="code" ,all.x=TRUE)
        sum(prev_count_hab$Freq * prev_count_hab$incl_prob)
    } else {
        NA
    }
}

hexa_area$prob_hab_qc <- pb_hab_hexas_qc
hexa_area$prob_hab_ca <- pb_hab_hexas_ca
# st_write(hexa_area, "outputs/schefferVille/hexas_shefferville.shp")
