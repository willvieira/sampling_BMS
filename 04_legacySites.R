########################################################################
### Legacy sites
### Author: Will Vieira
### February 28, 2021
########################################################################

########!IMPORTANT! ########
# Update on 28 jan 2023
# This script is no longer used as now all legacy sites are provided by the user in a single csv file on script `08_runGRTS.R`
############################

library(tidyverse)
library(sf)


########################################################################
# Steps
#   - Load spatial data and legacy site points
#   - Filter legacy sites points that can be used
#   - Count number of points by hexagon
########################################################################


# load all info
districts <- readRDS('data/districts.RDS')
legacy_sites <- sf::st_read('rawLayers/Registre_LegacyV1_W4.shp') %>%
                sf::st_transform(
                    sf::st_crs(
                        sf::st_read(
                            grep(
                                'phab',
                                list.files(
                                    'output',
                                    recursive = TRUE,
                                    pattern = '.shp',
                                    full.names = TRUE
                                ),
                                value = TRUE
                            )[1]
                        )
                    )
                )

# Filter points for BMS
legacy_sites <- subset(legacy_sites, BMS == 'OUI')


# Count legacy_sites per hexagon for each ecoregion
for(id in unique(districts$ECOREGION))
{
    # load hexagons
    hexa_ecoregion <- sf::st_read(paste0('output/', tolower(gsub(' ', '_', unique(subset(districts, ECOREGION == id)$REGION_NAM))), '_', id, '/hexa_cost.shp'), quiet = TRUE)
    
    # count number of lagacy site points
    hexa_ecoregion$legacySite <- lengths(sf::st_intersects(hexa_ecoregion, legacy_sites))

    # save as new shapefile
    sf::write_sf(hexa_ecoregion,
                 paste0('output/',
                        tolower(gsub(' ', '_', unique(subset(districts, ECOREGION == id)$REGION_NAM))),
                        '_', id, '/hexa_', id, '.shp'))

}
