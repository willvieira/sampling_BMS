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
doParallel::registerDoParallel(6)

load("data/spatialVectors.rda")
land_ca <- raster("data/landcover_ca_30m.tif")
#land_qc <- raster("data/landcover_qc_30m.tif")



###############################################
# Habitat prevalence within each ecoregions ---------------------------------------------------------
###############################################
# See page 32: Protocols/BMS_Boreal Monit_LB_DRAFT_May2_2018.docx

# Function to compute the prevalence habitat (area by habitat / area district) for each ecodistricts
prevalence_hab <- function(district, landcover)
{
    # Subselect area
    landcover <- mask(crop(landcover, district), district)

    # Inclusion prob
    # Count pixels by habitat type
    freq_hab <- as.data.frame(table(landcover[]))
    names(freq_hab) <- c("code", "freq")

    # Compute inclusion probability
    freq_hab$incl_prob <- 1 / (nrow(freq_hab)) / freq_hab$freq

    # Add poly ID
    freq_hab$ID_poly <- district$ECODISTRIC

    # Return the dataframe
    return(freq_hab)
}


# Canada
prev_ca_by_district <- list()
for(i in seq_len(nrow(districts)))
{
    prev_ca_by_district[[i]] <- prevalence_hab(district = districts[i,], landcover = land_ca)
    cat('   Computing prevalence ', round(i/nrow(districts) * 100, 2), '%\r')
}

# # Quebec
# prev_qc_by_district <- list()
# for(i in seq_len(nrow(districts)))
# {
#     prev_qc_by_district[[i]] <- prevalence_hab(district = districts[i,], landcover = land_ca)
#     print(i)
# }



####### Save prevalence by district
prev_all_ca <- do.call(rbind, prev_ca_by_district[which(sapply(prev_ca_by_district, class) == "data.frame")])
saveRDS(prev_all_ca, 'data/prev_all_ca.RDS')
#prev_all_qc <- do.call(rbind, prev_qc_by_district[which(sapply(prev_qc_by_district, class) == "data.frame")])
#saveRDS(prev_all_qc, 'data/prev_all_qc.RDS')


# viz how probability vary in function of # of veg type (code) and its frequence
par(mfrow = c(1, 2))
# qc
nbCode <- seq(11, 49)
freqCode <- seq(min(prev_all_qc$freq), quantile(prev_all_qc$freq, .75), by = 2000)
z <- outer(nbCode, freqCode, function(nbCode, freqCode) 1/nbCode/freqCode)
fields::image.plot(nbCode, freqCode, z, xlab = '# Vegetation type', ylab = 'Veg Freq')
points(prev_all_qc$code, prev_all_qc$freq, cex = 0.5, pch = 19, col = rgb(0, 0, 0, 0.6))
# ca
nbCode <- seq(6, 11)
freqCode <- seq(min(prev_all_ca$freq), quantile(prev_all_ca$freq, .75), by = 2000)
z <- outer(nbCode, freqCode, function(nbCode, freqCode) 1/nbCode/freqCode)
fields::image.plot(nbCode, freqCode, z, xlab = 'Vegetation type', ylab = '')
points(prev_all_ca$code, prev_all_ca$freq, cex = 0.5, pch = 19, col = rgb(0, 0, 0, 0.6))





###############################################
# Select ecoregion ---------------------------------------------------------
###############################################

prev_all_ca <- readRDS('data/prev_all_ca.RDS')

for (id in districts$ECODISTRIC) {

    print(glue('District {id} -- Process {which(districts$ECODISTRIC == id)} on {length(districts$ECODISTRIC)}'))
    # Select district
    district <- districts[districts$ECODISTRIC == id,]

    # Select hexagons within the area (st_insersects preserve topology)
    hexa_district <- hexa[which(st_intersects(hexa, district, sparse = FALSE)), ]
    
    # Select prevalence for that district
    #prev_qc_district <- subset(prev_all_qc, ID_poly == id)
    prev_ca_district <- subset(prev_all_ca, ID_poly == id)

    # Extract habitat values within each hexagons -- FOR QC
    #print(glue('Start to extract habitat values for QC'))

    #hab_qc <- foreach(i = seq_len(nrow(hexa_district)), .packages = "raster", .combine = c) %dopar% {    
    #    count_qc <- as.data.frame(table(raster::extract(land_qc, hexa_district[i,])))
    #    if(nrow(count_qc) > 0){
    #        prev_count_hab <- merge(count_qc, prev_qc_district, by.x = "Var1", by.y = "code" , all.x = TRUE)
    #        sum(prev_count_hab$Freq * prev_count_hab$incl_prob)
    #    } else {
    #        NA
    #    }
    #}

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
    #hexa_district$hab_qc <- hab_qc

    # Add new phab columns in attributes table
    print(glue('Saving'))

    # Writing results in corresponding folder
    path <- glue('outputs/{tolower(str_replace_all(unique(district$REGION_NAM), " ", "_"))}/{id}/')
    dir.create(path, recursive = TRUE, showWarnings = FALSE)

    write_sf(hexa_district, paste0(path,"hexa_phab.shp"))
}

stopImplicitCluster()
