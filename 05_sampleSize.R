########################################################################
### Calculate sample size for each ecoregion
### Author: Will Vieira
### February 9, 2021
########################################################################


library(tidyverse)
library(sf)


ecoregions <- readRDS('data/ecoregions.RDS')
hexas <- readRDS('data/hexa_complete.RDS')


########################################################################
# Steps
# - For each ecoregion
#   - Get hexagons in which centroid are inside ecoregion
#   - Union selected hexagons together into a single polygon
#   - Calculate total area of the polygon/ecoregion
# - Calculate sample size based on polygon/ecoregion area
# - save object with ecoregion name/code and sample size
########################################################################



out_dt <- data.frame()
for(eco in ecoregions)
{
    # get hexagons for specific region
    hexa_eco <- subset(hexas, ecoregion == eco)

    # filter NA
    hexa_eco <- hexa_eco[which(!is.na(hexa_eco$propNA)), ]
    
    # filter propNA <= 80%
    hexa_ecoProp <- subset(hexa_eco, propNA <= 0.8)

    # output
    out_dt <- rbind(out_dt, data.frame(ecoregion = eco, hexagonTotal = nrow(hexa_eco), hexagonProp80 = nrow(hexa_ecoProp), legacySites = sum(hexa_eco$legacySite), area = st_area(st_union(hexa_eco))))
}



# Calculate sample size based on polygon/ecoregion area (removing of ecoregions that have been divided)
eco_ignore <- which(out_dt$ecoregion %in% gsub('N', '', grep('N', out_dt$ecoregion, value = TRUE)))
totalArea <- sum(out_dt$area[-eco_ignore])
out_dt$sizeProp <- out_dt$area/totalArea


# save
saveRDS(out_dt, file = 'data/nbHexa_ecoregion.RDS')
