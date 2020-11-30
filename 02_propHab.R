########################################################################
### Title: Compute inclusion probability based on habitat prevalence
### Author: Steve Vissault
### Last edited by Will Vieira
########################################################################


library(raster)
library(sf)
library(glue)
library(stringr)
library(doParallel)
doParallel::registerDoParallel(parallel::detectCores() - 2)

load("data/spatialVectors.rda")
land_ca <- raster("data/landcover_ca_30m.tif")
#land_qc <- raster("data/landcover_qc_30m.tif")



###############################################
# Habitat prevalence within each ecoregions ---------------------------------------------------------
###############################################
# See page 32: Protocols/BMS_Boreal Monit_LB_DRAFT_May2_2018.docx

# Function to compute the prevalence habitat (area by habitat / area ecoregion) for each ecoregion
prevalence_hab <- function(ecoregion, landcover)
{
    # Subselect area
    landcover <- mask(crop(landcover, ecoregion), ecoregion)

    # Inclusion prob
    # Count pixels by habitat type
    freq_hab <- as.data.frame(table(landcover[]))
    names(freq_hab) <- c("code", "freq")

    # Compute inclusion probability
    freq_hab$incl_prob <- 1 / (nrow(freq_hab)) / freq_hab$freq

    # Add ecoregion ID
    freq_hab$ID_ecoregion <- rep(unique(ecoregion$ECOREGION), nrow(freq_hab))

    # Return the dataframe
    return(freq_hab)
}


# Using land Canada
###########################

# Same as above but in parallel
prev_ca_by_ecoregion <- foreach(i = unique(districts$ECOREGION), .packages = "raster") %dopar% {
    try(prevalence_hab(ecoregion = subset(districts, ECOREGION == i), landcover = land_ca))
}


# # Using land Quebec
###########################

# prev_qc_by_ecoregion <- foreach(i = unique(districts$ECOREGION), .packages = "raster") %dopar% {
#     try(prevalence_hab(ecoregion = subset(districts, ECOREGION == i), landcover = land_qc))
# }


####### Save prevalence by district
prev_all_ca <- do.call(rbind, prev_ca_by_ecoregion[which(sapply(prev_ca_by_ecoregion, class) == "data.frame")])
saveRDS(prev_all_ca, 'data/prev_all_ca.RDS')
#prev_all_qc <- do.call(rbind, prev_qc_by_ecoregion[which(sapply(prev_qc_by_ecoregion, class) == "data.frame")])
#saveRDS(prev_all_qc, 'data/prev_all_qc.RDS')


# # viz how probability vary in function of # of veg type (code) and its frequence
# par(mfrow = c(1, 2))
# # qc
# nbCode <- seq(11, 49)
# freqCode <- seq(min(prev_all_qc$freq), quantile(prev_all_qc$freq, .75), by = 2000)
# z <- outer(nbCode, freqCode, function(nbCode, freqCode) 1/nbCode/freqCode)
# fields::image.plot(nbCode, freqCode, z, xlab = '# Vegetation type', ylab = 'Veg Freq')
# points(prev_all_qc$code, prev_all_qc$freq, cex = 0.5, pch = 19, col = rgb(0, 0, 0, 0.6))
# # ca
# nbCode <- seq(6, 11)
# freqCode <- seq(min(prev_all_ca$freq), quantile(prev_all_ca$freq, .75), by = 2000)
# z <- outer(nbCode, freqCode, function(nbCode, freqCode) 1/nbCode/freqCode)
# fields::image.plot(nbCode, freqCode, z, xlab = 'Vegetation type', ylab = '')
# points(prev_all_ca$code, prev_all_ca$freq, cex = 0.5, pch = 19, col = rgb(0, 0, 0, 0.6))





###############################################
# Select ecoregion ---------------------------------------------------------
###############################################



prev_all_ca <- readRDS('data/prev_all_ca.RDS')

for (id in unique(districts$ECOREGION)) {

    print(glue('Ecoregion {id} -- Process {which(unique(districts$ECOREGION) == id)} on {length(unique(districts$ECOREGION))}'))
    
    # Select ecoregion
    ecoregion <- districts[districts$ECOREGION == id, ]

    # Select hexagons within the area (st_insersects preserve topology)
    hexa_ecoregion <- hexa[which(st_intersects(hexa, st_union(ecoregion), sparse = FALSE)), ]
    
    # Select prevalence for that district
    #prev_qc_ecoregion <- subset(prev_all_qc, ID_ecoregion == id)
    prev_ca_ecoregion <- subset(prev_all_ca, ID_ecoregion == id)

    # Extract habitat values within each hexagons -- FOR QC
    #print(glue('Start to extract habitat values for QC'))

    #hab_qc <- foreach(i = seq_len(nrow(hexa_ecoregion)), .packages = "raster", .combine = c) %dopar% {    
    #    count_qc <- as.data.frame(table(raster::extract(land_qc, hexa_ecoregion[i, ])))
    #    if(nrow(count_qc) > 0){
    #        prev_count_hab <- merge(count_qc, prev_qc_ecoregion, by.x = "Var1", by.y = "code" , all.x = TRUE)
    #        sum(prev_count_hab$Freq * prev_count_hab$incl_prob)
    #    } else {
    #        NA
    #    }
    #}

    # Extract habitat values within each hexagons -- FOR CA
    print(glue('Start to extract habitat values for CA'))
    hab_ca <- foreach(i = seq_len(nrow(hexa_ecoregion)), .packages = "raster", .combine = c) %dopar% {    
        count_ca <- as.data.frame(table(extract(land_ca, hexa_ecoregion[i, ])))
        if(nrow(count_ca) > 0)
        {
            if(sum(count_ca$Freq) > 1000) # ensure that hexagon has at least 100 pixels for secondary sampling
            {
            prev_count_hab <- merge(count_ca, prev_ca_ecoregion, by.x="Var1", by.y="code" ,all.x=TRUE)
            sum(prev_count_hab$Freq * prev_count_hab$incl_prob)
        } else {
            NA
        }
        } else {
            NA
        }
    }

    # Assigning new columns with hab prob
    hexa_ecoregion$hab_ca <- hab_ca
    #hexa_ecoregion$hab_qc <- hab_qc

    # Add new phab columns in attributes table
    print(glue('Saving'))

    # Writing results in corresponding folder
    path <- glue('output/{tolower(str_replace_all(unique(ecoregion$REGION_NAM), " ", "_"))}_{id}/')
    dir.create(path, recursive = TRUE, showWarnings = FALSE)

    write_sf(hexa_ecoregion, paste0(path, 'hexa_phab.shp'))
}

stopImplicitCluster()
