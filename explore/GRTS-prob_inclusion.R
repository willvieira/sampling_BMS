# Run GRTS for each layer of probability

library(raster)
library(exactextractr)
library(tidyverse)
library(spatstat)
library(sf)
library(spsurvey)
library(MBHdesign)

set.seed(0.0)



# parameters

    # Maximum proportion of NA pixels (non-habitat) in a hexagon
    # if 0.8, it means an hexagon has to have at least 20% of habitat pixels
    prop_na = 0.8

    # Percentage of hexagons to cover a region
    sample_effort = 0.02

    # Number of replications when running the GRTS
    nb_rep = 1

    # Ecoregions
    eco_sim = '101'
    
    # Name of file and columns to extract legacy
    legacyFile = file.path('data', 'legacySites.csv')
    lat = 'latitude'
    lon = 'longitude'

    # Buffer size (in Km) to adjust inclusion probability of hexagons around legacy sites
    bufferSize_p = 10

    # Buffer size (in Km) to adjust sample size of ecoregion in function of the number and distribution of legacy site
    # It can be a number that will be used for all ecoregion or a file address with the specific buffer size for each ecoregion
    bufferSize_N = file.path('data', 'bufferSize_N.csv')

    # Output folder to save the shapefiles with PSU and SSU
    outputFolder = file.path('..', '..', 'ownCloud', 'BMS_Bruno', 'selection2023')

#



# Prepare hexagons

    hexas <- readRDS(file.path('data', 'hexa_complete.RDS')) %>%
        filter(propNA <= prop_na) %>%
        filter(ecoregion %in% eco_sim) %>%
        mutate(
            p = (hab_prob * cost_prob) / sum(hab_prob * cost_prob)
        ) %>%
        filter(p != 0)

#



# legacy sites

    # function to transform Latitude & longitude legacy site points in a table
    # with the number of points per hexagon ID (ET_Index)
    import_legacySites <- function(File, lat_name, lon_name)
    {
        hx <- hexas %>%
            st_transform(4326)
        
        # read file
        lg <- read_csv(File, show_col_types = FALSE) %>%
            rename(
                lat = all_of(lat_name),
                lon = all_of(lon_name)
            ) %>%
            st_as_sf(
                coords = c('lon', 'lat'),
                crs = st_crs(hx)
            )

        # intersect
        nbLegacy <- hx %>%
            st_contains(lg, sparse = FALSE) %>%
            apply(1, sum)

        tibble(
            ET_Index = hx$ET_Index,
            nbLegacySites = nbLegacy
        ) %>%
        filter(nbLegacySites > 0)
    }


    # load and transform legacy sites (slow function)
    legacySites <- import_legacySites(
        File = legacyFile,
        lat_name = lat,
        lon_name = lon
    )

    # merge to hexagons
    hexas <- hexas %>%
        left_join(legacySites) %>%
        mutate(nbLegacySites = replace_na(nbLegacySites, 0))

#



# Calculate sample size given number of hexagons and legacy sites

    # define buffer size to adjust sample size
    if(is.character(bufferSize_N)) {
        if(file.exists(bufferSize_N)) {
            buffSizeN <- read_csv(bufferSize_N, show_col_types = FALSE) |>
                mutate(ecoregion = as.character(ecoregion))
        }else{
            stop(paste0('File "', bufferSize_N, '" does not exist. Please check if the name is correct.'))
        }
    }else if(is.numeric(bufferSize_N)){
        buffSizeN <- tibble(
            ecoregion = eco_sim,
            bufferSize = bufferSize_N
        )
    }else{
        stop('Type of `bufferSize_N` must be either a numeric or a character')
    }
    
    # function to get sample size for a specific ecoregion given:
    # number of hexagons, number of legacy sites, bufferSize, and sample effort
    get_sampleSize <- function(eco, hx, bf_N, sample_e)
    {
        # get the hexagons centroid for a ecoregion
        hexa_eco <- hx %>%
            filter(ecoregion == eco) %>%
            st_centroid()

        # create a buffer of size BufferSize_N around legacy sites
        hexa_legacy_bf <- hexa_eco %>%
            filter(nbLegacySites > 0) %>%
            st_buffer(subset(bf_N, ecoregion == eco)$bufferSize * 1000) %>%
            st_union()

        # Compute number of hexagons in which centroid is inside legacy buffer
        nbHexas_legacy <- hexa_eco %>%
            st_intersects(hexa_legacy_bf) %>%
            unlist() %>%
            sum()
        
        # Compute adjusted sample size
        adj_sampleSize <- round((nrow(hexa_eco) - nbHexas_legacy) * sample_e, 0)

        # Only if ecoregion is too small, assure to sample at least two hexagons for all ecoregions
        if((nrow(hexa_eco) * sample_e) < 2 & nbHexas_legacy < 1)
            adj_sampleSize = 2

        return(adj_sampleSize)
    }


    sampleSize <- map_dbl(
        setNames(eco_sim, paste0('eco_', eco_sim)),
        get_sampleSize,
        hx = hexas,
        bf_N = buffSizeN,
        sample_e = sample_effort
    )

#
 


# GRTS simulations

    # Sample size by stratum (ecoregion)
    Stratdsgn  <- sampleSize[sampleSize > 0]

    # make sure that ecoregion has at least 3 times more hexagons than sample N
    eco_to_remove <- hexas |>
        st_drop_geometry() |>
        group_by(ecoregion) |>
        summarise(nbHexas = n()) |>
        left_join(
            sampleSize |>
                enframe() |>
                mutate(
                    ecoregion = as.character(parse_number(name)),
                    sampleSize_extra = value * 3
                ) |>
                select(ecoregion, sampleSize_extra)
        ) |>
        mutate(
            diff = nbHexas - sampleSize_extra
        ) |>
        filter(diff < 0) |>
        pull(ecoregion)
    
    Stratdsgn <- Stratdsgn[!names(Stratdsgn) %in% paste0('eco_', eco_to_remove)]
        
    # Prepare sample frame
    sampleFrame <- hexas %>%
        filter(ecoregion %in% parse_number(names(Stratdsgn))) |>
        mutate(    
            eco_name = paste0('eco_', ecoregion), # to match design name
            mdcaty_hab  = sum(Stratdsgn) * hab_prob/sum(hab_prob),
            mdcaty_cost  = sum(Stratdsgn) * cost_prob/sum(cost_prob),
            mdcaty_p  = sum(Stratdsgn) * p/sum(p),
            geometry = sf::st_geometry(sf::st_centroid(geometry))
        )
        
    # Coordinates of all hexagons in matrix format for MBHdesign
    coord_mt <- sampleFrame %>%
        st_coordinates()

    # Coordinates of legacy hexagons
    legacySites <- sampleFrame %>%
        filter(nbLegacySites > 0) %>%
        st_coordinates()

    sampleFrame$adj_p <- MBHdesign::alterInclProbs(
        legacy.sites = legacySites,
        potential.sites = coord_mt,
        inclusion.probs = sampleFrame$mdcaty_p,
        sigma = bufferSize_p * 1000
    )

    # run GRTS
    ###############################################

    set.seed(10)
    
    # habitat
    out_sample_hab <- spsurvey::grts(
        sframe = sampleFrame,
        n_base = Stratdsgn,
        stratum_var = 'eco_name',
        aux_var = 'mdcaty_hab'
    )

    hexas$sel_hab <- 0
    hexas$sel_hab[which(hexas$ET_Index %in% out_sample_hab$sites_base$ET_Index)] <- 1

    # Cost
    out_sample_cost <- spsurvey::grts(
        sframe = sampleFrame,
        n_base = Stratdsgn,
        stratum_var = 'eco_name',
        aux_var = 'mdcaty_cost'
    )

    hexas$sel_cost <- 0
    hexas$sel_cost[which(hexas$ET_Index %in% out_sample_cost$sites_base$ET_Index)] <- 1

    # p
    out_sample_p <- spsurvey::grts(
        sframe = sampleFrame,
        n_base = Stratdsgn,
        stratum_var = 'eco_name',
        aux_var = 'mdcaty_p'
    )

    hexas$sel_p <- 0
    hexas$sel_p[which(hexas$ET_Index %in% out_sample_p$sites_base$ET_Index)] <- 1

    # legacy
    out_sample_leg <- spsurvey::grts(
        sframe = sampleFrame,
        n_base = sampleSize,
        stratum_var = 'eco_name',
        aux_var = 'adj_p'
    )

    hexas$sel_leg <- 0
    hexas$sel_leg[which(hexas$ET_Index %in% out_sample_leg$sites_base$ET_Index)] <- 1

#



# Viz

    p1 <- hexas |>
        filter(!is.na(hab_prob)) |>
        filter(hab_prob < quantile(hab_prob, probs = 0.999)) |>
        ggplot() +
            geom_sf(aes(fill = hab_prob), color = 'transparent') +
            theme_minimal() +
            scale_fill_viridis_c() +
            coord_sf(crs = 4326) +
            xlab('') + ylab('') +
            theme(legend.position = 'none') +
            labs(subtitle = 'Habitat')
    p2 <- hexas |>
        filter(!is.na(cost_prob)) |>
        ggplot() +
            geom_sf(aes(fill = cost_prob), color = 'transparent') +
            theme_minimal() +
            scale_fill_viridis_c() +
            coord_sf(crs = 4326) +
            xlab('') + ylab('') +
            theme(legend.position = 'none') +
            labs(subtitle = 'Coût')
    p3 <- hexas |>
        filter(!is.na(p)) |>
        filter(p < quantile(p, probs = 0.999)) |>
        ggplot() +
            geom_sf(aes(fill = p), color = 'transparent') +
            theme_minimal() +
            scale_fill_viridis_c() +
            coord_sf(crs = 4326) +
            xlab('') + ylab('') +
            theme(legend.position = 'none') +
            labs(subtitle = 'Habitat + coût')


    p4 <- hexas |>
        ggplot() +
            geom_sf(aes(fill = sel_hab), color = 'transparent') +
            theme_minimal() +
            scale_fill_viridis_c() +
            coord_sf(crs = 4326) +
            xlab('') + ylab('') +
            theme(legend.position = 'none') +
            labs(subtitle = '')
    p5 <- hexas |>
        ggplot() +
            geom_sf(aes(fill = sel_cost), color = 'transparent') +
            theme_minimal() +
            scale_fill_viridis_c() +
            coord_sf(crs = 4326) +
            xlab('') + ylab('') +
            theme(legend.position = 'none') +
            labs(subtitle = '')
    p6 <- hexas |>
        ggplot() +
            geom_sf(aes(fill = sel_p), color = 'transparent') +
            theme_minimal() +
            scale_fill_viridis_c() +
            coord_sf(crs = 4326) +
            xlab('') + ylab('') +
            theme(legend.position = 'none') +
            labs(subtitle = '')

    ggpubr::ggarrange(p1, p2, p3, p4, p5, p6, nrow = 2, ncol = 3)


    p7 <- hexas |>
        mutate(adj_p = sampleFrame$adj_p) |>
        filter(!is.na(adj_p)) |>
        filter(adj_p < quantile(adj_p, probs = 0.999)) |>
        ggplot() +
            geom_sf(aes(fill = adj_p), color = 'transparent') +
            #geom_sf(data = hexas |> filter(nbLegacySites > 0), fill = 'red', color = 'transparent') +
            theme_minimal() +
            scale_fill_viridis_c() +
            coord_sf(crs = 4326) +
            xlab('') + ylab('') +
            theme(legend.position = 'none') +
            labs(subtitle = 'Habitat + coût + sites historiques')

        
    p8 <- hexas |>
        ggplot() +
            geom_sf(aes(fill = sel_leg), color = 'transparent') +
            geom_sf(data = hexas |> filter(nbLegacySites > 0), fill = 'darkblue', color = 'transparent') +
            theme_minimal() +
            scale_fill_viridis_c() +
            coord_sf(crs = 4326) +
            xlab('') + ylab('') +
            theme(legend.position = 'none') +
            labs(subtitle = '')

    ggpubr::ggarrange(p3, p7, p6, p8, nrow = 2, ncol = 2)

#
