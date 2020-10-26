### Title: Compute inclusion probability based on habitat prevalence
### Author: Steve Vissault

library(raster)
library(sf)
library(glue)
library(stringr)

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
registerDoParallel(cores = detectCores()-2)
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

load("data/prevalence_by_district.rda")

for (id in districts$ECODISTRIC) {

    print(glue('District {id} -- Process {which(districts$ECODISTRIC == id)} on {length(districts$ECODISTRIC)}'))
    # Select district
    district <- districts[districts$ECODISTRIC == id,]

    # Select hexagons within the area (st_insersects preserve topology)
    hexa_district <- hexa[which(st_intersects(hexa, district, sparse = FALSE)),]
    
    # Select prevalence for that district
    prev_qc_district <- subset(prev_all_qc, ID_poly == id)
    prev_ca_district <- subset(prev_all_ca, ID_poly == id)

    # Extract habitat values within each hexagons -- FOR QC
    print(glue('Start to extract habitat values for QC'))

    hab_qc <- foreach(i = seq_len(nrow(hexa_district)), .packages = "raster", .combine = c) %dopar% {    
        count_qc <- as.data.frame(table(extract(land_qc, hexa_district[i,])))
        if(nrow(count_qc) > 0){
            prev_count_hab <- merge(count_qc, prev_qc_district, by.x="Var1", by.y="code" ,all.x=TRUE)
            sum(prev_count_hab$Freq * prev_count_hab$incl_prob)
        } else {
            NA
        }
    }

    # Extract habitat values within each hexagons -- FOR CA
    print(glue('Start to extract habitat values for CA'))
    hab_ca <- foreach(i = seq_len(nrow(hexa_district)), .packages = "raster", .combine = c) %dopar% {    
        count_ca <- as.data.frame(table(extract(land_ca, hexa_district[i,])))
        if(nrow(count_ca) > 0){
            prev_count_hab <- merge(count_ca, prev_ca_district, by.x="Var1", by.y="code" ,all.x=TRUE)
            sum(prev_count_hab$Freq * prev_count_hab$incl_prob)
        } else {
            NA
        }
    }

    # Assigning new columns with hab prob
    hexa_district$hab_ca <- hab_ca
    hexa_district$hab_qc <- hab_qc

    # Add new phab columns in attributes table
    print(glue('Saving'))

    # Writing results in corresponding folder
    path <- glue('outputs/{tolower(str_replace_all(unique(district$REGION_NAM), " ", "_"))}/{id}/')
    dir.create(path, recursive = TRUE, showWarnings = FALSE)

    write_sf(hexa_district, paste0(path,"hexa_phab.shp"))
}


