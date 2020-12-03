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
# Get sum of habitat probability for each hexagon using habitat prevalence from the ecoregion
###############################################

# load data calculated above
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
    hab_ca <- foreach(i = seq_len(nrow(hexa_ecoregion)), .packages = "raster") %dopar% {    
        values_ca <- extract(land_ca, hexa_ecoregion[i, ])[[1]]
        count_ca <- as.data.frame(table(values_ca))
        if(nrow(count_ca) > 0)
        {
            # save proportion of empty pixels in the hexagonq
            propNA <- sum(is.na(values_ca))/length(values_ca)
            # probability of inclusion
            prev_count_hab <- merge(count_ca, prev_ca_ecoregion, by.x = "values_ca", by.y = "code" , all.x = TRUE)
            # return
            c(propNA, sum(prev_count_hab$Freq * prev_count_hab$incl_prob))
        } else {
            c(NA, NA)
        }
    }

    # Assigning new columns with hab prob and proportion of NA in the hexagon attributes table
    hexa_ecoregion$hab_ca <- unlist(lapply(hab_ca, `[[`, 2))
    hexa_ecoregion$propNA <- unlist(lapply(hab_ca, `[[`, 1))

    print(glue('Saving'))

    # Writing results in corresponding folder
    path <- glue('output/{tolower(str_replace_all(unique(ecoregion$REGION_NAM), " ", "_"))}_{id}/')
    dir.create(path, recursive = TRUE, showWarnings = FALSE)

    write_sf(hexa_ecoregion, paste0(path, 'hexa_phab.shp'))
}

stopImplicitCluster()




###############################################
# Plot percentage of NA for each hexagon
###############################################

pdf('proportionNA_ecoregion.pdf', width = 12, height = 6)
for(id in unique(districts$ECOREGION))
{
    # Select ecoregion
    ecoregion <- districts[districts$ECOREGION == id, ]

    # load ecoregion hexagons
    hexa_ecoregion <- sf::st_read(paste0('output/', tolower(gsub(' ', '_', unique(ecoregion$REGION_NAM))), '_', id, '/hexa_phab.shp'), quiet = TRUE)

    # plot
    par(mfrow = c(1, 2), mar = c(2.5, 2.5, 2.5, 3.5), mgp = c(1.4, 0.2, 0), tck = -.008)
    cuts <- cut(seq(0, 100, 5), c(-Inf,50,Inf))
    hist(hexa_ecoregion$propNA * 100, breaks = seq(0, 100, 5), main = '', ylab = 'Frequency (# of hexagons)', xlab = 'Proportion of NA pixels within each hexagon (%)', col = c(rgb(100, 100, 100, 200, maxColorValue = 255), rgb(71, 141, 255, 200, maxColorValue = 255))[cuts])
    box()
    mtext(paste('Ecoregion', id, '-', unique(ecoregion$REGION_NAM)), 3, line = -2, outer = TRUE)
    
    hexa_ecoFilt <- subset(hexa_ecoregion, propNA >= 0.5)
    xLim <- c(0, max(hexa_ecoregion$hab_ca, na.rm = TRUE))
    h <- hist(hexa_ecoregion$hab_ca, main = '', ylab = 'Frequency (# of hexagons)', xlab = 'Habitat probability of hexagons', breaks = 30, col = rgb(100, 100, 100, 210, maxColorValue = 255))
    par(new = TRUE)
    hist(hexa_ecoFilt$hab_ca, breaks = h$breaks, xlim = xLim, xlab = '', ylab = '', main = '', xaxt = 'n', yaxt = 'n', col = rgb(116, 169, 255, 140, maxColorValue = 255))
    axis(4, col = rgb(71, 141, 255, maxColorValue = 255), col.axis = rgb(71, 141, 255, maxColorValue = 255))
    mtext('Frequency (hexagons with more than 50% of NA pixels', 4, line = 1, cex = 0.9, col = rgb(71, 141, 255, maxColorValue = 255))
    box()
}
dev.off()
