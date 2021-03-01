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
doParallel::registerDoParallel(parallel::detectCores() - 1)

load("data/spatialVectors.rda")
land_ca <- raster("data/landcover_ca_30m.tif")
land_caqc <- raster("data/landcover_caqc_30m.tif")



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

# Calculate habitat prevalence for each ecoregion in parallel
# Using two different habitat layers depending on the ecoregion

# Ecoregions to use land_ca (using land_ca because there is no information from land_caqc available)
eco_land_ca <- c(28, 30, 46, 48, 49, 78, 86, 216)
eco_land_caqc <- setdiff(unique(districts$ECOREGION), eco_land_ca)

prev_ca_by_ecoregion <- foreach(i = eco_land_ca, .packages = "raster") %dopar% {
    try(prevalence_hab(ecoregion = subset(districts, ECOREGION == i), landcover = land_ca))
}

# Adjust crs to the match with land_caqc
districts_caqc <- sf::st_transform(districts, sf::st_crs(land_caqc))

prev_caqc_by_ecoregion <- foreach(i = eco_land_caqc, .packages = "raster") %dopar% {
    try(prevalence_hab(ecoregion = subset(districts_caqc, ECOREGION == i), landcover = land_caqc))
}


####### Save prevalence by district
prev_all <- do.call(rbind, c(
                    prev_ca_by_ecoregion[which(sapply(prev_ca_by_ecoregion, class) == "data.frame")],
                    prev_caqc_by_ecoregion[which(sapply(prev_caqc_by_ecoregion, class) == "data.frame")]))

saveRDS(prev_all, 'data/prev_all.RDS')


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
# freqCode <- seq(min(prev_all$freq), quantile(prev_all$freq, .75), by = 2000)
# z <- outer(nbCode, freqCode, function(nbCode, freqCode) 1/nbCode/freqCode)
# fields::image.plot(nbCode, freqCode, z, xlab = 'Vegetation type', ylab = '')
# points(prev_all$code, prev_all$freq, cex = 0.5, pch = 19, col = rgb(0, 0, 0, 0.6))





###############################################
# Get sum of habitat probability for each hexagon using habitat prevalence from the ecoregion
###############################################


# Create object with hexagon centroid (it will be used to filter hexagons within an ecoregion)
hexa_cent <- hexa
hexa_cent$geometry <- hexa_cent %>%
                        sf::st_centroid() %>%
                        sf::st_geometry() %>%
                        sf::st_transform(st_crs(districts))

# Temporary hexa with crs ajusted to match land_caqc
hexa_caqc <- hexa %>%
                sf::st_transform(crs(land_caqc))

# Get habitat codes
habCodes <- sort(as.numeric(levels(prev_all$code)))

for (id in unique(districts$ECOREGION))
{
    print(glue('Ecoregion {id} -- Process {which(unique(districts$ECOREGION) == id)} on {length(unique(districts$ECOREGION))}'))
    
    # Select ecoregion
    ecoregion <- districts[districts$ECOREGION == id, ]

    # Select prevalence for that district
    #prev_qc_ecoregion <- subset(prev_all_qc, ID_ecoregion == id)
    prev_all_ecoregion <- subset(prev_all, ID_ecoregion == id)

    # Select hexagons within the ecoregion (st_insersects preserve topology)
    hexaPos <- which(st_intersects(hexa_cent, st_union(ecoregion), sparse = FALSE))
    hexa_ecoregion <- hexa[hexaPos, ]
    hexa_ecoregion_caqc <- hexa_caqc[hexaPos, ]

    # Extract habitat values within each hexagons
    print(glue('Start to extract habitat values for CA'))

    if(id %in% eco_land_ca)
    {        
        hab_ls <- foreach(i = seq_len(nrow(hexa_ecoregion)), .packages = "raster") %dopar% {    
            values_hab <- raster::extract(land_ca, hexa_ecoregion[i, ])[[1]]
            count_hab <- as.data.frame(table(values_hab))
            if(nrow(count_hab) > 0)
            {
                # save proportion of empty pixels in the hexagonq
                propNA <- sum(is.na(values_hab))/length(values_hab)
                # probability of inclusion
                prev_count_hab <- merge(count_hab, prev_all_ecoregion, by.x = "values_hab", by.y = "code" , all.x = TRUE)
                # habitat frequency
                habFreq <- setNames(rep(0, length(habCodes)), habCodes)
                habFreq[as.character(count_hab$values_hab)] <- count_hab$Freq

                # return
                c(propNA = propNA,
                hab_prob = sum(prev_count_hab$Freq * prev_count_hab$incl_prob),
                setNames(habFreq, paste0('land_ca_', names(habFreq))))

            } else {
                rep(NA, 23)
            }
        }
    }else if(id %in% eco_land_caqc)
    {
        hab_ls <- foreach(i = seq_len(nrow(hexa_ecoregion_caqc)), .packages = "raster") %dopar% {    
            values_hab <- raster::extract(land_caqc, hexa_ecoregion_caqc[i, ])[[1]]
            count_hab <- as.data.frame(table(values_hab))
            if(nrow(count_hab) > 0)
            {
                # save proportion of empty pixels in the hexagonq
                propNA <- sum(is.na(values_hab))/length(values_hab)
                # probability of inclusion
                prev_count_hab <- merge(count_hab, prev_all_ecoregion, by.x = "values_hab", by.y = "code" , all.x = TRUE)
                # habitat frequency
                habFreq <- setNames(rep(0, length(habCodes)), habCodes)
                habFreq[as.character(count_hab$values_hab)] <- count_hab$Freq

                # return
                c(propNA = propNA,
                hab_prob = sum(prev_count_hab$Freq * prev_count_hab$incl_prob),
                setNames(habFreq, paste0('land_ca_', names(habFreq))))

            } else {
                rep(NA, 23)
            }
        }
    }else{
        warning('Ecoregion ', id, ' was not found in either `eco_land_ca` nor `eco_land_caqc. Skipping...`')
    }

    # Assign new columns with (i) hab prob, (ii) proportion of NA, and (iii) frequency of habitats
    # in the hexagon attributes table
    hab_dt <- do.call(rbind, hab_ls)
    hexa_ecoregion <- cbind(hexa_ecoregion, hab_dt)

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
    h <- hist(hexa_ecoregion$propNA * 100, breaks = seq(0, 100, 5), main = '', ylab = 'Frequency (# of hexagons)', xlab = 'Proportion of NA pixels within each hexagon (%)', col = c(rgb(100, 100, 100, 200, maxColorValue = 255), rgb(71, 141, 255, 200, maxColorValue = 255))[cuts])
    legend('topright', legend = paste0(round(sum(h$counts[h$mids > 50])/sum(h$counts)*100, 1), '% of hexagons have a\nproportion > 50% of NA pixels'), cex = .9, bty = 'n')
    box()
    mtext(paste('Ecoregion', id, '-', unique(ecoregion$REGION_NAM)), 3, line = -2, outer = TRUE)
    
    hexa_ecoFilt <- subset(hexa_ecoregion, propNA >= 0.5)
    xLim <- c(0, max(hexa_ecoregion$hab_prob, na.rm = TRUE))
    h <- hist(hexa_ecoregion$hab_prob, main = '', ylab = 'Frequency (# of hexagons)', xlab = 'Habitat probability of hexagons', breaks = 30, col = rgb(100, 100, 100, 210, maxColorValue = 255))
    par(new = TRUE)
    hist(hexa_ecoFilt$hab_prob, breaks = h$breaks, xlim = xLim, xlab = '', ylab = '', main = '', xaxt = 'n', yaxt = 'n', col = rgb(116, 169, 255, 140, maxColorValue = 255))
    axis(4, col = rgb(71, 141, 255, maxColorValue = 255), col.axis = rgb(71, 141, 255, maxColorValue = 255))
    mtext('Frequency (hexagons with more than 50% of NA pixels', 4, line = 1, cex = 0.9, col = rgb(71, 141, 255, maxColorValue = 255))
    box()

    # par(las = 2)
    # land_dt <- st_drop_geometry(hexa_ecoregion[, paste0('land_ca_', c(1:14, 16))])
    # landTotal <- apply(land_dt, 2, sum, na.rm = TRUE)
    # landPropNA <- apply(land_dt[which(hexa_ecoregion$propNA >= 0.5), ], 2, sum, na.rm = TRUE)
    # barplot(landTotal, horiz = TRUE)
    # barplot(landPropNA, add = TRUE, col = rgb(71, 141, 255, maxColorValue = 255), horiz = TRUE)
    # paste0(round(landPropNA/landTotal * 100, 1), '%')

}
dev.off()
