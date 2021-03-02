########################################################################
### Calculate sample size for each ecoregion
### Author: Will Vieira
### February 9, 2021
########################################################################


library(tidyverse)
library(sf)

load('data/spatialVectors.rda')


########################################################################
# Steps
# - For each ecoregion
#   - Get hexagons in which centroid are inside ecoregion
#   - Union selected hexagons together into a single polygon
#   - Calculate total area of the polygon/ecoregion
# - Calculate sample size based on polygon/ecoregion area
# - save object with ecoregion name/code and sample size
########################################################################


# load hexagons
hexa_ls <- list()
for (id in unique(districts$ECOREGION))
{
  hexa_ls[[as.character(id)]] <- sf::st_read(paste0('output/', tolower(gsub(' ', '_', unique(subset(districts, ECOREGION == id)$REGION_NAM))), '_', id, '/hexa_', id, '.shp'), quiet = TRUE)
  hexa_ls[[as.character(id)]]$ecoregion <- id
}

# rbind all shapefiles into one sf oject
hexas <- do.call(rbind, hexa_ls)
rm(hexa_ls)


# Hexagon centroid
hexa_cen <- hexas %>%
            st_centroid() %>%
            st_transform(st_crs(districts))


# unique ecoregion names and codes
ecoregions <- setNames(districts$ECOREGION, districts$REGION_NAM)


out_dt <- data.frame()
for(eco in ecoregions)
{
    # ecoregion polygon
    eco_p <- subset(districts, ECOREGION == eco)
    
    # get hexagons for specific region
    hexa_eco <- hexas[which(st_intersects(hexa_cen, st_union(eco_p), sparse = FALSE)), ]

    # filter NA
    hexa_eco <- hexa_eco[which(!is.na(hexa_eco$propNA)), ]
    
    # filter propNA <= 80%
    hexa_ecoProp <- subset(hexa_eco, propNA <= 0.8)

    # output
    out_dt <- rbind(out_dt, data.frame(ecoregion = eco, hexagonTotal = nrow(hexa_eco), hexagonProp80 = nrow(hexa_ecoProp), legacySites = sum(hexa_eco$legacySite), area = st_area(st_union(hexa_eco))))
}



# Calculate sample size based on polygon/ecoregion area
out_dt$sizeProp <- out_dt$area/sum(out_dt$area)


# save
saveRDS(out_dt, file = 'data/nbHexa_ecoregion.RDS')
