########################################################################
### Title: Run GRTS
### Author: Will Vieira
### December 9, 2020
### Last edited: Feb 5, 2023
########################################################################

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

    # Total sample size for SSU (Main + Over)
    # It must be an even number
    ssu_N = 6

    # Number of replications when running the GRTS
    nb_rep = 15

    # Ecoregions
    eco_sim = c(
        '7', '28', '30', '31',
        '46', '47', '48', '49',
        '73', '77', '78', '86',
        '72', '74', '75', '76',
        '96', #'96N', '96S',
        '99',
        '100', #'100N', '100S',
        '101', #'101N', '101S',
        '102',
        '103', #'103N', '103S',
        '117',
        '216',
        '217'#, '217N', '217S'
    )
    
    # Name of file and columns to extract legacy
    legacyFile = file.path('data', 'legacySites.csv')
    lat = 'latitude'
    lon = 'longitude'

    # Buffer size (in Km) to adjust inclusion probability of hexagons around legacy sites
    bufferSize_p = 10

    # Buffer size (in Km) to adjust sample size of ecoregion in function of the number and distribution of legacy site
    # It can be a number that will be used for all ecoregion or a file address with the specific buffer size for each ecoregion
    bufferSize_N = file.path('data', 'bufferSize_N.csv')

    # Distance between SSU centroid (in meters)
    ssu_dist = 294


    # Output folder to save the shapefiles with PSU and SSU
    outputFolder = file.path('..', '..', 'ownCloud', 'BMS_Bruno', 'selection2023')

    # suffix to add for each output layer
    # e.g.: PSU-SOQB_ALL-SUFFIX.shp
    fileSuffix = 'V2023'

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

    # Prepare sample frame
    sampleFrame <- hexas %>%
        mutate(    
            eco_name = paste0('eco_', ecoregion), # to match design name
            mdcaty  = sum(Stratdsgn) * p/sum(p),
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
        inclusion.probs = sampleFrame$mdcaty,
        sigma = bufferSize_p * 1000
    )

    # run GRTS    
    grts_out <- list()
    for(Rep in 1:nb_rep)
    {
        out_sample <- spsurvey::grts(
            sframe = sampleFrame,
            n_base = Stratdsgn,
            stratum_var = 'eco_name',
            aux_var = 'adj_p'
        )

        grts_out[[paste0('Rep_', Rep)]] <- out_sample$sites_base$ET_Index
    }

#



# Select cheapest replication

    cheapest_rep <- map_df(
        grts_out,
        ~ hexas %>%
            st_drop_geometry() %>%
            filter(ET_Index %in% .x) %>%
            summarise(totalCost = sum(costSum))
    ) %>%
    pull(totalCost) %>%
    which.min()


    selected_hexas <- hexas %>%
        filter(ET_Index %in% grts_out[[cheapest_rep]])

#



## EXTRA SAMPLE
# Select the closest hexagon with the highest inclusion probability

    # loop over ecoregions to make sure neighbours are from the same ecoregion
    selected_extra <- hexas[0, ]
    for(eco in eco_sim)
    {
        hexas_eco <- subset(hexas, ecoregion == eco)
        hexas_sel_eco <- subset(selected_hexas, ecoregion == eco)

        # Extract neighbours hexagons
        neigh_mt <- hexas_eco %>%
            st_centroid() %>%
            st_intersects(
                y = st_buffer(hexas_sel_eco, dist = 3500),
                sparse = FALSE
            )

        # Remove the focus hexagon (the selected one)
        rr = match(hexas_sel_eco$ET_Index, hexas_eco$ET_Index)
        cc = seq_along(rr)
        neigh_mt[rr + nrow(neigh_mt) * (cc - 1)] <- FALSE

        # Select the extra hexagon based on the highest p
        best_p <- apply(
            neigh_mt,
            2,
            function(x)
                hexas_eco$ET_Index[x][which.max(hexas_eco$p[x])]
        )

        # if a selected hexagon has no neighbour, select a random from the ecoregion
        toCheck <- unlist(lapply(best_p, length))
        if(any(toCheck == 0)) {
            # which hexagons were not selected? 
            nonSelected_hexas <- setdiff(hexas_eco$ET_Index, hexas_sel_eco$ET_Index)
            
            # sample from non selected hexas
            best_p[which(toCheck == 0)] <- sample(nonSelected_hexas, sum(toCheck == 0))
        }
        
        selected_extra <- rbind(
            selected_extra,
            hexas_eco[match(best_p, hexas_eco$ET_Index), ]
        )
    }

#



# Select SSU

    # Calculate probability following habitat proportion from ecoregion
    ########################

    land_ca <- raster('data/landcover_ca_30m.tif')
    prev_all <- readRDS('data/prev_all.RDS')

    # function to generate SSU points (from: https://github.com/dhope/BASSr)
    genSSU <- function(h, spacing)
    {
        ch <- as_tibble(st_coordinates(h))
        top_point <- ch[which.max(ch$Y),]
        bottom_point <- ch[which.min(ch$Y),]
        gridsize <- 2*floor(abs(top_point$Y-bottom_point$Y)/spacing)+3
        rowAngle <- tanh((top_point$X-bottom_point$X)/(top_point$Y-bottom_point$Y))


        cent <- st_centroid(h) %>%
            bind_cols(as_tibble(st_coordinates(.))) %>%
            st_drop_geometry %>%
            dplyr::select(ET_Index, X, Y)


        genRow <- function(cX, cY, sp,...){
            tibble(rowid = seq(-gridsize,gridsize)) %>%
            mutate(X = sin(60*pi/180+rowAngle) *sp*rowid + {{cX}},
                    Y = cos(60*pi/180+rowAngle) *sp*rowid  + {{cY}})
        }

        centroids <- tibble(crowid=seq(-gridsize,gridsize)) %>%
            mutate(cY = cos(rowAngle) *spacing*crowid + cent$Y,
                #spacing/2*crowid + cent$Y,
                cX =  sin(rowAngle) *spacing*crowid + cent$X) %>%
            #cent$X + crowid* sqrt(spacing**2-(spacing/2)**2)) %>%
            rowwise() %>%
            mutate(row = list(genRow(cX = cX,cY = cY,sp = spacing))) %>%
            unnest(row) %>%
            dplyr::select(X,Y) %>%
            st_as_sf(coords = c("X", "Y"), crs = st_crs(h)) %>%
            st_filter(h) %>%
            mutate(
                ET_Index = h$ET_Index,
                ecoregion = h$ecoregion,
                ssuID = row_number()
            )
        return(centroids)
    }


    # function to sample SSU
    sample_SSU <- function(ssuid, prob, geom, filtered, ssuDist, N)
    {
        # check if N is even
        if(N %% 2 != 0)
            stop('`ssu_N` must be a even number.')

        filtered_1 <- filtered
        out <- rep(0, length(filtered))

        # loop to sample 4
        count = 1
        while(
            count <= N &
            sum(get(paste0('filtered_', count))) > 1
        ){
            # sample point
            assign(
                paste0('sample_', count),
                sample(
                    ssuid[get(paste0('filtered_', count))],
                    size = 1,
                    prob = prob[get(paste0('filtered_', count))]
                )
            )

            # remove points around the first sample for second point
            toKeep <- !st_intersects(
                geom,
                st_buffer(
                    geom[which(get(paste0('sample_', count)) == ssuid)],
                    dist = ssuDist * 2 + ssuDist * 0.1),
                    sparse = FALSE
            )[, 1]

            # update available points
            assign(
                paste0('filtered_', count + 1),
                get(paste0('filtered_', count)) & toKeep
            )

            # assign point to output vector
            out[get(paste0('sample_', count))] = count
            
            count = count + 1
        }

        return( out )

    }



    # Generate SSU points
    SSUs <- map_dfr(
        seq_len(nrow(selected_hexas)),
        ~ genSSU(selected_hexas[.x, ], spacing = ssu_dist)
    )

    # Buffer of half `ssu_dist` to compute habitat inclusion prob
    SSU_bf <- st_buffer(SSUs, dist = ssu_dist/2)

    # extract pixels for each SSU polygon
    hab_pixels <- exactextractr::exact_extract(
        land_ca,
        SSU_bf,
        progress = FALSE
    )
    rm(SSU_bf)

    # get frequence of each class of habitat
    count_hab <- Map(
                function(x, y) {
                    freq <- table(x$value)
                    if(length(freq) > 0) {
                        data.frame(freq, ecoregion = y)
                    }else{
                        NA
                    }
                },
                x = hab_pixels,
                y = SSUs$ecoregion
            )

    # merge with inclusion probability
    # and calculate inclusion probbaility for each NON empty polygon
    SSUs$incl_prob <- unlist(
            lapply(
                count_hab,
                function(x) {
                    if(is.data.frame(x)) {
                        mg_df <- merge(x, subset(prev_all, ID_ecoregion == x$ecoregion[1]), by.x = "Var1", by.y = "code" , all.x = TRUE)
                        sum(mg_df$Freq * mg_df$incl_prob)
                    }else{
                        NA
                    }
                }
            )
        )
    rm(count_hab)

    SSUs <- SSUs %>%
        group_by(ET_Index) %>%
        mutate(
            incl_prob = incl_prob/sum(incl_prob, na.rm = TRUE)
        ) %>%
        ungroup()


    # Calculate proportion of NA
    SSUs$propNA <- map_dbl(
        hab_pixels,
        ~ sum(is.na(.x$value))/nrow(.x)
    )
    rm(hab_pixels)

    # neighbours matrix
    neighbour_ls <- list()
    for(hx in unique(SSUs$ET_Index))
    {
        ssuhx <- subset(SSUs, ET_Index == hx)

        neighbour_ls[[hx]] <- st_intersects(
            ssuhx,
            st_buffer(ssuhx, dist = ssu_dist + ssu_dist * 0.1),
            sparse = FALSE
        )
    }


    # These are the following rules to a SSU be available to be sampled:
    # - Must have 6 neighbours (less than that means it's a border SSU)
    # - Must have at least 1 - `prop_na` of non empty pixels
    # - 4 out 6 neighbours must also respect the above rule
    SSU_selected = SSUs %>%
        group_by(ET_Index) %>%
        mutate(
            nbNeighb = map_dbl(
                ssuID,
                ~ sum(neighbour_ls[[unique(ET_Index)]][, .x]) - 1
            ),
            nbPropNA = map_int(
                ssuID,
                .f = function(x, pNA, ETI)
                    sum(
                        pNA[setdiff(which(neighbour_ls[[ETI]][, x]), x)] <= prop_na
                    ),
                pNA = propNA,
                ETI = unique(ET_Index)
            ),
            sampled = sample_SSU(
                ssuid = ssuID,
                prob = incl_prob,
                geom = geometry,
                filtered = propNA <= prop_na & nbNeighb == 6  & nbPropNA >= 4,
                ssuDist = ssu_dist,
                N = ssu_N
            )
        ) %>%
        ungroup()

    # Prepare selected SSU and their specific neighbours with code like `A_B`
    # `A` is for the SSU sample ID (1:`ssu_N`)
    # `B` is for the neighbour ID (0:6 where zero is the focal point, and 1:6 are the 6 neighbours starting from the top point moving clockwise)
    SSUmain <- subset(SSU_selected, sampled > 0)
    SSU_main_ls <- list()

    for(i in 1:nrow(SSUmain))
    {
        # get neighbours for specific row
        nei_hx <- SSU_selected %>% 
            filter(ET_Index == SSUmain$ET_Index[i]) %>% 
            filter(ssuID %in% which(neighbour_ls[[SSUmain$ET_Index[i]]][, SSUmain$ssuID[i]]))

        # code A 
        nei_hx$codeA <- SSUmain$sampled[i]
        
        # code B
        nei_hx$codeB <- c(4, 3, 5, 0, 2, 6, 1)

        SSU_main_ls[[i]] <- nei_hx
    }

    SSUmain <- do.call(rbind, SSU_main_ls)




    # SAME BUT FOR THE EXTRA HEXAGONS (TODO)
    # Generate SSU points
    SSUs <- map_dfr(
        seq_len(nrow(selected_extra)),
        ~ genSSU(selected_extra[.x, ], spacing = ssu_dist)
    )

    # Buffer of half `ssu_dist` to compute habitat inclusion prob
    SSU_bf <- st_buffer(SSUs, dist = ssu_dist/2)

    # extract pixels for each SSU polygon
    hab_pixels <- exactextractr::exact_extract(
        land_ca,
        SSU_bf,
        progress = FALSE
    )
    rm(SSU_bf)

    # get frequence of each class of habitat
    count_hab <- Map(
                function(x, y) {
                    freq <- table(x$value)
                    if(length(freq) > 0) {
                        data.frame(freq, ecoregion = y)
                    }else{
                        NA
                    }
                },
                x = hab_pixels,
                y = SSUs$ecoregion
            )

    # merge with inclusion probability
    # and calculate inclusion probbaility for each NON empty polygon
    SSUs$incl_prob <- unlist(
            lapply(
                count_hab,
                function(x) {
                    if(is.data.frame(x)) {
                        mg_df <- merge(x, subset(prev_all, ID_ecoregion == x$ecoregion[1]), by.x = "Var1", by.y = "code" , all.x = TRUE)
                        sum(mg_df$Freq * mg_df$incl_prob)
                    }else{
                        NA
                    }
                }
            )
        )
    rm(count_hab)

    SSUs <- SSUs %>%
        group_by(ET_Index) %>%
        mutate(
            incl_prob = incl_prob/sum(incl_prob, na.rm = TRUE)
        ) %>%
        ungroup()


    # Calculate proportion of NA
    SSUs$propNA <- map_dbl(
        hab_pixels,
        ~ sum(is.na(.x$value))/nrow(.x)
    )
    rm(hab_pixels)

    # neighbours matrix
    neighbour_ls <- list()
    for(hx in unique(SSUs$ET_Index))
    {
        ssuhx <- subset(SSUs, ET_Index == hx)

        neighbour_ls[[hx]] <- st_intersects(
            ssuhx,
            st_buffer(ssuhx, dist = ssu_dist + ssu_dist * 0.1),
            sparse = FALSE
        )
    }


    # These are the following rules to a SSU be available to be sampled:
    # - Must have 6 neighbours (less than that means it's a border SSU)
    # - Must have at least 1 - `prop_na` of non empty pixels
    # - 4 out 6 neighbours must also respect the above rule
    SSU_selected_extra = SSUs %>%
        group_by(ET_Index) %>%
        mutate(
            nbNeighb = map_dbl(
                ssuID,
                ~ sum(neighbour_ls[[unique(ET_Index)]][, .x]) - 1
            ),
            nbPropNA = map_int(
                ssuID,
                .f = function(x, pNA, ETI)
                    sum(
                        pNA[setdiff(which(neighbour_ls[[ETI]][, x]), x)] <= prop_na
                    ),
                pNA = propNA,
                ETI = unique(ET_Index)
            ),
            sampled = sample_SSU(
                ssuid = ssuID,
                prob = incl_prob,
                geom = geometry,
                filtered = propNA <= prop_na & nbNeighb == 6  & nbPropNA >= 4,
                ssuDist = ssu_dist,
                N = ssu_N
            )
        ) %>%
        ungroup()
    

    # Prepare selected SSU and their specific neighbours with code like `A_B`
    # `A` is for the SSU sample ID (1:`ssu_N`)
    # `B` is for the neighbour ID (0:6 where zero is the focal point, and 1:6 are the 6 neighbours starting from the top point moving clockwise)
    SSUover <- subset(SSU_selected_extra, sampled > 0)
    SSU_over_ls <- list()

    for(i in 1:nrow(SSUover))
    {
        # get neighbours for specific row
        nei_hx <- SSU_selected_extra %>% 
            filter(ET_Index == SSUover$ET_Index[i]) %>% 
            filter(ssuID %in% which(neighbour_ls[[SSUover$ET_Index[i]]][, SSUover$ssuID[i]]))

        # code A 
        nei_hx$codeA <- SSUover$sampled[i]
        
        # code B
        nei_hx$codeB <- c(4, 3, 5, 0, 2, 6, 1)

        SSU_over_ls[[i]] <- nei_hx
    }

    SSUover <- do.call(rbind, SSU_over_ls)



    # Prepare export of shapefiles by ecoregion
    for(eco in eco_sim)
    {
        # create folder
        eco_path <- file.path(outputFolder, paste0('ecoregion_', eco))
        dir.create(eco_path, recursive = TRUE)

        # PSU
        ###########################################

        varsToRm = c('OBJECTID', 'Join_Count', 'TARGET_FID', 'ET_ID', 'ET_ID_Old', 'ET_IDX_Old')

        hexas_eco <- hexas %>%
            filter(ecoregion == eco) %>%
            select(-all_of(varsToRm))

        # rename attributes table so it has a maximum of 10 characters
        names(hexas_eco) <- abbreviate(names(hexas_eco), minlength = 10)

        coords <- hexas_eco %>%
            sf::st_centroid() %>%
            sf::st_transform(4326) %>%
            sf::st_coordinates() %>%
            as.data.frame()

        # all hexagons
        hexas_eco %>%
            mutate(
                latitude = coords$Y,
                longitude = coords$X
            ) %>%
            write_sf(
                file.path(
                    eco_path,
                    paste0('PSU-SOBQ_eco', eco, '_ALL-', fileSuffix, '.shp')
                )
            )

        # main hexagons
        hexas_eco %>%
            mutate(
                latitude = coords$Y,
                longitude = coords$X
            ) %>%
            filter(ET_Index %in% subset(selected_hexas, ecoregion == eco)$ET_Index) %>%
            write_sf(
                file.path(
                    eco_path,
                    paste0('PSU-SOBQ_eco', eco, '_Main-', fileSuffix, '.shp')
                )
            )

        # over hexagons
        hexas_eco %>%
            mutate(
                latitude = coords$Y,
                longitude = coords$X
            ) %>%
            filter(ET_Index %in% subset(selected_extra, ecoregion == eco)$ET_Index) %>%
            write_sf(
                file.path(
                    eco_path,
                    paste0('PSU-SOBQ_eco', eco, '_Over-', fileSuffix, '.shp')
                )
            )

        
        # SSU
        ###########################################

        # SSU main
        SSU_eco <- SSU_selected %>%
            filter(ecoregion == eco) %>%
            select(-c('nbNeighb', 'nbPropNA', 'sampled'))

        coords <- SSU_eco %>%
            sf::st_transform(4326) %>%
            sf::st_coordinates() %>%
            as.data.frame()

        SSU_eco %>%
            mutate(
                latitude = coords$Y,
                longitude = coords$X
            ) %>%
            write_sf(
                file.path(
                    eco_path,
                    paste0('SSU-SOBQ_eco', eco, '_ALL_main-', fileSuffix, '.shp')
                )
            )

        # SSU over
        SSU_eco_extra <- SSU_selected_extra %>%
            filter(ecoregion == eco) %>%
            select(-c('nbNeighb', 'nbPropNA', 'sampled'))

        coords <- SSU_eco_extra %>%
            sf::st_transform(4326) %>%
            sf::st_coordinates() %>%
            as.data.frame()

        SSU_eco_extra %>%
            mutate(
                latitude = coords$Y,
                longitude = coords$X
            ) %>%
            write_sf(
                file.path(
                    eco_path,
                    paste0('SSU-SOBQ_eco', eco, '_ALL_over-', fileSuffix, '.shp')
                )
            )

        # SSU selected main
        SSUmain_eco <- SSUmain %>%
            filter(ecoregion == eco) %>%
            mutate(
                sample = paste0(codeA, '_', codeB)
            ) %>%
            select('ET_Index', 'ecoregion', 'ssuID', 'incl_prob', 'sample')

        coords <- SSUmain_eco %>%
            sf::st_transform(4326) %>%
            sf::st_coordinates() %>%
            as.data.frame()
        
        SSUmain_eco %>%
            mutate(
                latitude = coords$Y,
                longitude = coords$X
            ) %>%
            write_sf(
                file.path(
                    eco_path,
                    paste0('SSU-SOBQ_eco', eco, '_selected_main-', fileSuffix, '.shp')
                )
            )
        
        # SSU selected over
        SSUover_eco <- SSUover %>%
            filter(ecoregion == eco) %>%
            mutate(
                sample = paste0(codeA, '_', codeB)
            ) %>%
            select('ET_Index', 'ecoregion', 'ssuID', 'incl_prob', 'sample')

        coords <- SSUover_eco %>%
            sf::st_transform(4326) %>%
            sf::st_coordinates() %>%
            as.data.frame()
        
        SSUover_eco %>%
            mutate(
                latitude = coords$Y,
                longitude = coords$X
            ) %>%
            write_sf(
                file.path(
                    eco_path,
                    paste0('SSU-SOBQ_eco', eco, '_selected_over-', fileSuffix, '.shp')
                )
            )
    }
