########################################################################
### Split ecoregions crossing the 50.3" latitude into North and South
### Due to bias from legacy sites that were limited to the 50.3"
### Author: Will Vieira
### March 31, 2021
########################################################################


library(tidyverse)
library(sf)

load('data/spatialVectors.rda')



########################################################################
# Steps
# - For each ecoregion
#   - Get centroid of hexagons
#   - Split hexagons whithin ecoregion between above and bellow 50.3"
#   - Save new hexagons as a new ecoregion ecoX_S and ecoX_S
# - Add the new ecoregions to the list of ecoregions to be used later
########################################################################



ecoregions_to_split <- c(96, 100, 101, 103, 217)


for(eco in ecoregions_to_split)
{
    # load hexagons from ecoregion
    hexa_ecoregion <- sf::read_sf(paste0('output/',
                                        tolower(gsub(' ', '_', unique(subset(districts, ECOREGION == eco)$REGION_NAM))),
                                        '_', eco, '/hexa_', eco, '.shp'))

     # Split ecoregion in two parts
    coords <- hexa_ecoregion %>%
                sf::st_centroid() %>%
                 sf::st_transform(4326) %>%
                sf::st_coordinates() %>%
                as.data.frame()

    hexa_ecoregionN <- hexa_ecoregion[coords$Y >= 50.3, ]
    hexa_ecoregionS <- hexa_ecoregion[!(coords$Y >= 50.3), ]
    
    # save new splited ecoregion
    sapply(paste0('output/', tolower(gsub(' ', '_', unique(subset(districts, ECOREGION == eco)$REGION_NAM))), '_', eco, c('N', 'S')), dir.create)

    write_sf(hexa_ecoregionN,
             paste0('output/', tolower(gsub(' ', '_', unique(subset(districts, ECOREGION == eco)$REGION_NAM))), '_', eco, 'N/hexa_', eco, 'N.shp'))
    write_sf(hexa_ecoregionS,
             paste0('output/', tolower(gsub(' ', '_', unique(subset(districts, ECOREGION == eco)$REGION_NAM))), '_', eco, 'S/hexa_', eco, 'S.shp'))             
}


# Create list of ecoregion to be used later
ecoregions <- sort(setNames(unique(districts$ECOREGION), tolower(gsub(' ', '_', unique(districts$REGION_NAM)))))
names(ecoregions_to_split) <- lapply(ecoregions_to_split, function(x) tolower(gsub(' ', '_', unique(subset(districts, ECOREGION == x)$REGION_NAM))))
ecoregions <- c(ecoregions,
                setNames(paste0(rep(ecoregions_to_split, each = 2), c('N', 'S')),
                         rep(names(ecoregions_to_split), each = 2)))

# sort list
eco_sort <- gsub('N', '.1', ecoregions)
eco_sort <- gsub('S', '.2', eco_sort)
ecoregions <- ecoregions[match(sort(as.numeric(eco_sort)), as.numeric(eco_sort))]

# Save with existing objs
save(area, districts, hexa, roads, trails, aeroports, train, ecoregions, file = "data/spatialVectors.rda")
